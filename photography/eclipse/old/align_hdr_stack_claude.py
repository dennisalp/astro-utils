#!/usr/bin/env python3
"""
align_hdr_stack.py
===================
Align and blend a bracketed sequence of solar-eclipse exposures (16-bit
TIFFs exported from Lightroom) into a single, well-exposed image.

Written for sequences like:
    7 frames, 1.7 EV apart, 0.5 s ... 1/2000 s
    Sony a7C II + Tamron 150-500 mm

------------------------------------------------------------------------
PIPELINE
------------------------------------------------------------------------
1. Load all TIFFs, kept in 16-bit for maximum quality.
2. Align every frame to a reference frame (default: the middle exposure)
   using ECC (Enhanced Correlation Coefficient) registration on an image
   pyramid. Frames are NOT aligned directly to the reference if they are
   several stops away from it -- instead alignment is *chained* through
   neighbouring exposures (1.7 EV apart, so lots of shared structure)
   and the resulting transforms are composed. This is far more robust
   than trying to match, say, the 1/2000 s frame directly against the
   0.5 s frame, which may barely share any visible structure.
3. Merge the aligned frames using one of:
     - "mertens" (default): exposure fusion. No exposure times or
       camera response curve needed, very robust, and the direct output
       is already a displayable image. This is effectively the same
       algorithm as enfuse/Hugin and is the standard approach for solar
       corona/prominence bracket stacks.
     - "debevec": "true" radiometric HDR merge using your exposure
       times, producing a linear-light HDR radiance map that is then
       tone-mapped back down to a displayable image. More faithful to
       actual relative brightness, but more fiddly, and OpenCV's
       implementation requires converting to 8-bit internally.
4. Write a 16-bit TIFF (and optionally an 8-bit preview), carrying over
   the ICC colour profile of the first input file if one is present.

------------------------------------------------------------------------
IMPORTANT CAVEATS -- PLEASE READ
------------------------------------------------------------------------
* Exposure times: for --method debevec you must supply the *actual*
  shutter speeds used, in seconds, in the same order as the (sorted)
  input files. A helper (--auto-exposures) can generate a theoretical
  geometric sequence from a base time and an EV step, but your camera
  almost certainly rounded to its nearest 1/3-stop shutter speeds, so
  prefer reading the real values from the EXIF of your original ARW
  files over trusting the theoretical sequence.
* Lightroom TIFF exports are usually NOT scene-linear (a tone curve /
  camera profile has been applied), which technically violates the
  reciprocity assumption behind Debevec-style HDR merging. In practice
  OpenCV's CalibrateDebevec estimates a response curve empirically from
  your own bracket, which absorbs most of this non-linearity, but for
  eclipse work "mertens" is the more predictable and widely-used choice
  and is the default here.
* Alignment model: solar-eclipse brackets shot in a few seconds on a
  tripod/star-tracker usually only need translation. --align euclidean
  additionally recovers rotation (e.g. field rotation on an alt-az
  mount) at the cost of being slightly less robust to converge.
* Memory: a 33 MP camera (a7C II) at 7 frames easily needs a few GB of
  RAM in 16-bit. Use --max-dimension to downscale for a quick preview
  run first if you're short on memory.

------------------------------------------------------------------------
DEPENDENCIES
------------------------------------------------------------------------
    pip install opencv-python numpy tifffile

------------------------------------------------------------------------
USAGE
------------------------------------------------------------------------
    # Simplest: point at a folder of 7 TIFFs, Mertens fusion, translation align
    python align_hdr_stack.py -i ./tiffs -o eclipse_fused.tif

    # Euclidean alignment (handles slight rotation) + explicit reference frame
    python align_hdr_stack.py -i ./tiffs -o eclipse_fused.tif \\
        --align euclidean --ref-index 3

    # "True" HDR merge with your real shutter speeds, Reinhard tonemap
    python align_hdr_stack.py -i ./tiffs -o eclipse_hdr.tif \\
        --method debevec \\
        --exposures 0.5 0.16667 0.05 0.0166667 0.005 0.0015625 0.0005 \\
        --tonemap reinhard

    # Also write an 8-bit JPEG preview alongside the 16-bit TIFF
    python align_hdr_stack.py -i ./tiffs -o eclipse_fused.tif --preview
"""

import argparse
import glob
import os
import sys

import cv2
import numpy as np
import tifffile

ICC_TAG_CODE = 34675  # standard TIFF tag for an embedded ICC profile

MOTION_TYPES = {
    "translation": cv2.MOTION_TRANSLATION,
    "euclidean": cv2.MOTION_EUCLIDEAN,
    "affine": cv2.MOTION_AFFINE,
}


# --------------------------------------------------------------------------
# I/O helpers
# --------------------------------------------------------------------------

def find_input_files(input_arg):
    """Accept either a directory (scanned for tif/tiff) or an explicit glob."""
    if os.path.isdir(input_arg):
        files = sorted(
            glob.glob(os.path.join(input_arg, "*.tif"))
            + glob.glob(os.path.join(input_arg, "*.tiff"))
            + glob.glob(os.path.join(input_arg, "*.TIF"))
            + glob.glob(os.path.join(input_arg, "*.TIFF"))
        )
    else:
        files = sorted(glob.glob(input_arg))
    # de-duplicate while preserving order (case-insensitive filesystems can
    # otherwise return the same file via two of the patterns above)
    seen = set()
    unique_files = []
    for f in files:
        rp = os.path.realpath(f)
        if rp not in seen:
            seen.add(rp)
            unique_files.append(f)
    return unique_files


def read_icc_profile(path):
    """Best-effort extraction of an embedded ICC profile from a TIFF."""
    try:
        with tifffile.TiffFile(path) as tf:
            tags = tf.pages[0].tags
            if ICC_TAG_CODE in tags:
                return tags[ICC_TAG_CODE].value
    except Exception as e:
        print(f"  [info] could not read ICC profile from {path}: {e}")
    return None


def load_images(files, max_dimension=None):
    """Load TIFFs as a list of HxWx3 uint16 arrays (RGB order)."""
    images = []
    shape = None
    for f in files:
        arr = tifffile.imread(f)
        if arr.ndim == 2:
            arr = np.stack([arr] * 3, axis=-1)
        if arr.shape[-1] == 4:
            arr = arr[..., :3]  # drop alpha if present
        if arr.dtype != np.uint16:
            # Handle 8-bit source gracefully by promoting to 16-bit range
            if arr.dtype == np.uint8:
                arr = (arr.astype(np.uint16)) * 257
            else:
                raise ValueError(
                    f"{f}: unexpected dtype {arr.dtype}, expected uint16 (or uint8)"
                )
        if max_dimension is not None:
            h, w = arr.shape[:2]
            scale = max_dimension / max(h, w)
            if scale < 1.0:
                arr = cv2.resize(
                    arr, (int(round(w * scale)), int(round(h * scale))),
                    interpolation=cv2.INTER_AREA,
                )
        if shape is None:
            shape = arr.shape
        elif arr.shape != shape:
            raise ValueError(
                f"{f}: shape {arr.shape} does not match first image shape {shape}. "
                "All frames must be the same resolution/orientation."
            )
        images.append(arr)
        print(f"  loaded {f}  ({arr.shape[1]}x{arr.shape[0]}, {arr.dtype})")
    return images


def save_tiff_16bit(path, image_float01, icc_profile=None):
    """image_float01: HxWx3 float in [0,1] -> written as 16-bit TIFF."""
    out = np.clip(image_float01, 0.0, 1.0)
    out_u16 = (out * 65535.0 + 0.5).astype(np.uint16)
    extratags = []
    if icc_profile is not None:
        extratags.append((ICC_TAG_CODE, "B", len(icc_profile), icc_profile, True))
    tifffile.imwrite(path, out_u16, photometric="rgb", extratags=extratags)


def save_preview(path, image_float01):
    out = np.clip(image_float01, 0.0, 1.0)
    out_u8 = (out * 255.0 + 0.5).astype(np.uint8)
    # cv2 expects BGR for imwrite; our arrays are RGB throughout
    cv2.imwrite(path, cv2.cvtColor(out_u8, cv2.COLOR_RGB2BGR))


# --------------------------------------------------------------------------
# Alignment
# --------------------------------------------------------------------------

def to_gray_f32(img_u16):
    """Rec.709 luma, normalized to roughly [0,1] float32 for ECC."""
    f = img_u16.astype(np.float32)
    gray = 0.2126 * f[..., 0] + 0.7152 * f[..., 1] + 0.0722 * f[..., 2]
    gray /= 65535.0
    return gray


def build_pyramid(img, levels):
    pyr = [img]
    for _ in range(levels - 1):
        pyr.append(cv2.pyrDown(pyr[-1]))
    return pyr[::-1]  # coarsest first


def ecc_align_pair(template_gray, input_gray, motion_type, levels=4,
                    iterations=300, eps=1e-7):
    """
    Returns the 2x3 warp matrix M such that:
        cv2.warpAffine(input_image, M, size,
                        flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)
    produces `input_image` resampled into `template`'s coordinate frame.

    Falls back to the identity transform (with a warning) if ECC fails to
    converge at any pyramid level, rather than crashing the whole run.
    """
    template_pyr = build_pyramid(template_gray, levels)
    input_pyr = build_pyramid(input_gray, levels)
    warp = np.eye(2, 3, dtype=np.float32)
    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, iterations, eps)
    for lvl in range(levels):
        t = template_pyr[lvl]
        i = input_pyr[lvl]
        try:
            _, warp = cv2.findTransformECC(t, i, warp, motion_type, criteria, None, 5)
        except cv2.error as e:
            print(f"    [warn] ECC did not converge at pyramid level {lvl}: {e}")
            print(f"    [warn] falling back to best warp found so far for this pair")
            break
        if lvl < levels - 1:
            warp[:, 2] *= 2.0  # translation scales when moving to a finer level
    return warp


def _to_3x3(m2x3):
    m = np.eye(3, dtype=np.float64)
    m[:2, :] = m2x3
    return m


def _to_2x3(m3x3):
    return m3x3[:2, :].astype(np.float32)


def align_stack(images, ref_index, motion_type, levels=4):
    """
    Align every frame in `images` to `images[ref_index]`.

    Alignment is chained through neighbouring exposures (adjacent bracket
    steps share the most common, well-exposed detail) and the pairwise
    transforms are composed to map every frame directly into the
    reference frame's coordinate system. Returns a list of 2x3 warp
    matrices (one per frame; identity for the reference itself), ready
    to be used with cv2.warpAffine(..., flags=INTER_LINEAR+WARP_INVERSE_MAP).
    """
    n = len(images)
    grays = [to_gray_f32(im) for im in images]

    cumulative = [None] * n
    cumulative[ref_index] = np.eye(3, dtype=np.float64)

    # frames below the reference: chain k -> k+1 -> ... -> ref
    for k in range(ref_index - 1, -1, -1):
        print(f"  aligning frame {k} -> {k + 1}")
        pairwise = ecc_align_pair(grays[k + 1], grays[k], motion_type, levels)
        cumulative[k] = _to_3x3(pairwise) @ cumulative[k + 1]

    # frames above the reference: chain k -> k-1 -> ... -> ref
    for k in range(ref_index + 1, n):
        print(f"  aligning frame {k} -> {k - 1}")
        pairwise = ecc_align_pair(grays[k - 1], grays[k], motion_type, levels)
        cumulative[k] = _to_3x3(pairwise) @ cumulative[k - 1]

    return [_to_2x3(c) for c in cumulative]


def apply_warps(images, warps):
    aligned = []
    for img, warp in zip(images, warps):
        h, w = img.shape[:2]
        out = cv2.warpAffine(
            img, warp, (w, h),
            flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP,
            borderMode=cv2.BORDER_REPLICATE,
        )
        aligned.append(out)
    return aligned


# --------------------------------------------------------------------------
# Merging
# --------------------------------------------------------------------------

def merge_mertens(images_u16):
    """
    Exposure fusion. NOTE: OpenCV's MergeMertens internally treats pixel
    values as if they were on an 8-bit (0..255) scale regardless of the
    input array's actual dtype -- it does NOT auto-detect 16-bit range.
    Feeding it images normalized to [0,1] (as you might naively do for
    16-bit data) silently produces a near-black, broken result. We
    therefore rescale our 16-bit data onto a 0..255-equivalent float32
    range before calling it.
    """
    imgs_255 = [im.astype(np.float32) * (255.0 / 65535.0) for im in images_u16]
    merger = cv2.createMergeMertens()
    fused = merger.process(imgs_255)  # float32, nominally in [0,1]
    return np.clip(fused, 0.0, 1.0)


def merge_debevec(images_u16, exposure_times_s, tonemap_method="reinhard", gamma=1.0):
    """
    "True" radiometric HDR merge. OpenCV's CalibrateDebevec/MergeDebevec
    hard-require 8-bit (CV_8U) input, so we necessarily lose precision
    within each frame here (the whole point of the technique is
    recovering dynamic range ACROSS frames, so this is standard and
    expected). Returns a tone-mapped float32 image in [0,1].
    """
    imgs_8u = [
        np.clip(im.astype(np.float32) / 257.0, 0, 255).astype(np.uint8)
        for im in images_u16
    ]
    times = np.array(exposure_times_s, dtype=np.float32)

    calibrate = cv2.createCalibrateDebevec()
    response = calibrate.process(imgs_8u, times)

    merger = cv2.createMergeDebevec()
    hdr = merger.process(imgs_8u, times, response)  # float32 radiance, unbounded

    if tonemap_method == "reinhard":
        tm = cv2.createTonemapReinhard(gamma=gamma)
    elif tonemap_method == "drago":
        tm = cv2.createTonemapDrago(gamma=gamma)
    elif tonemap_method == "mantiuk":
        tm = cv2.createTonemapMantiuk(gamma=gamma)
    else:
        tm = cv2.createTonemap(gamma=gamma)

    ldr = tm.process(hdr)
    return np.clip(ldr, 0.0, 1.0)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def auto_exposure_sequence(n, base_seconds, ev_step):
    return [base_seconds * (2.0 ** (-ev_step * i)) for i in range(n)]


def parse_args():
    p = argparse.ArgumentParser(
        description="Align and HDR-stack a bracketed solar eclipse TIFF sequence.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True,
                    help="Directory of TIFFs, or a glob pattern (e.g. './tiffs/*.tif')")
    p.add_argument("-o", "--output", required=True,
                    help="Output 16-bit TIFF path")
    p.add_argument("--method", choices=["mertens", "debevec"], default="mertens",
                    help="Merge method (default: mertens)")
    p.add_argument("--align", choices=["translation", "euclidean", "affine", "none"],
                    default="translation", help="Alignment motion model (default: translation)")
    p.add_argument("--ref-index", type=int, default=None,
                    help="Index (0-based, after sorting) of the reference frame. "
                         "Default: the middle frame.")
    p.add_argument("--pyramid-levels", type=int, default=4,
                    help="Number of pyramid levels for ECC alignment (default: 4)")
    p.add_argument("--max-dimension", type=int, default=None,
                    help="Optionally downscale so the longer side is at most this "
                         "many pixels (useful for a fast preview run)")
    p.add_argument("--exposures", type=float, nargs="+", default=None,
                    help="Exposure times in seconds, in file order (required for "
                         "--method debevec). Example: 0.5 0.153 0.047 0.0146 "
                         "0.00449 0.00138 0.00043")
    p.add_argument("--auto-exposures", action="store_true",
                    help="Auto-generate a theoretical exposure sequence instead of "
                         "--exposures, using --base-exposure and --ev-step. Prefer "
                         "real EXIF values when possible.")
    p.add_argument("--base-exposure", type=float, default=0.5,
                    help="Longest exposure in seconds, for --auto-exposures (default: 0.5)")
    p.add_argument("--ev-step", type=float, default=1.7,
                    help="EV step between brackets, for --auto-exposures (default: 1.7)")
    p.add_argument("--tonemap", choices=["reinhard", "drago", "mantiuk", "linear"],
                    default="reinhard", help="Tonemap operator for --method debevec")
    p.add_argument("--preview", action="store_true",
                    help="Also write an 8-bit .jpg preview next to the output")
    return p.parse_args()


def main():
    args = parse_args()

    files = find_input_files(args.input)
    if len(files) < 2:
        sys.exit(f"Found {len(files)} TIFF file(s) at '{args.input}' -- need at least 2.")
    print(f"Found {len(files)} input files:")
    for f in files:
        print(f"  {f}")

    icc_profile = read_icc_profile(files[0])
    if icc_profile is not None:
        print("Embedded ICC profile found in first file; will carry it over to the output.")

    print("\nLoading images...")
    images = load_images(files, max_dimension=args.max_dimension)
    n = len(images)

    ref_index = args.ref_index if args.ref_index is not None else n // 2
    if not (0 <= ref_index < n):
        sys.exit(f"--ref-index {ref_index} out of range for {n} images")
    print(f"\nReference frame: index {ref_index} ({files[ref_index]})")

    # Validate/resolve exposure times up front (before the expensive alignment
    # step) so a --exposures mistake fails fast.
    exposures = None
    if args.method == "debevec":
        if args.auto_exposures:
            exposures = auto_exposure_sequence(n, args.base_exposure, args.ev_step)
            print("  using auto-generated exposure sequence (seconds):",
                  [f"{e:.5f}" for e in exposures])
        elif args.exposures is not None:
            if len(args.exposures) != n:
                sys.exit(f"--exposures has {len(args.exposures)} values but there are "
                          f"{n} input images.")
            exposures = args.exposures
        else:
            sys.exit("--method debevec requires --exposures (or --auto-exposures).")

    if args.align == "none":
        print("\nSkipping alignment (--align none).")
        aligned = images
    else:
        print(f"\nAligning ({args.align})...")
        motion_type = MOTION_TYPES[args.align]
        warps = align_stack(images, ref_index, motion_type, levels=args.pyramid_levels)
        aligned = apply_warps(images, warps)

    print(f"\nMerging ({args.method})...")
    if args.method == "mertens":
        fused = merge_mertens(aligned)
    else:
        fused = merge_debevec(aligned, exposures, tonemap_method=args.tonemap)

    print(f"\nSaving {args.output} ...")
    save_tiff_16bit(args.output, fused, icc_profile=icc_profile)

    if args.preview:
        preview_path = os.path.splitext(args.output)[0] + "_preview.jpg"
        save_preview(preview_path, fused)
        print(f"Saved preview {preview_path}")

    print("Done.")


if __name__ == "__main__":
    main()
