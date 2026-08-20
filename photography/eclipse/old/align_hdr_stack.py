#!/usr/bin/env python3
"""
align_hdr_stack.py
===================
Align and blend a bracketed sequence of solar-eclipse exposures (16-bit
TIFFs exported from Lightroom) into a single, well-exposed image.

------------------------------------------------------------------------
PIPELINE
------------------------------------------------------------------------
1. Load all TIFFs and automatically extract EXIF metadata (exposure
   times, camera model, lens, ISO, aperture).
2. Align every frame to a reference frame using ECC chained registration.
3. Merge the aligned frames using either Exposure Fusion (Mertens) or 
   radiometric HDR (Debevec).
4. Write a 16-bit TIFF carrying over the ICC colour profile and any 
   shared camera EXIF metadata (camera, lens, ISO, etc.).

------------------------------------------------------------------------
DEPENDENCIES
------------------------------------------------------------------------
    pip install opencv-python numpy tifffile exifread

------------------------------------------------------------------------
USAGE
------------------------------------------------------------------------
    # Simple Exposure Fusion (reads everything automatically)
    python align_hdr_stack.py -i ./tiffs -o eclipse_fused.tif

    # "True" HDR merge (reads actual exposure times from EXIF)
    python align_hdr_stack.py -i ./tiffs -o eclipse_hdr.tif \
        --method debevec --tonemap reinhard

    # Optional: Also write an 8-bit JPEG preview alongside the 16-bit TIFF
    python align_hdr_stack.py -i ./tiffs -o eclipse_fused.tif --preview
"""

import argparse
import glob
import os
import sys

import cv2
import numpy as np
import tifffile
import exifread

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


def get_file_metadata(filepath):
    """Extracts relevant EXIF data to be used in processing and output."""
    with open(filepath, 'rb') as f:
        tags = exifread.process_file(f, details=False)
        
    meta = {}
    
    # 1. Exposure Time (Crucial for Debevec)
    for tag in ['EXIF ExposureTime', 'Image ExposureTime']:
        if tag in tags:
            val = tags[tag].values[0]
            if hasattr(val, 'num') and hasattr(val, 'den') and val.den != 0:
                meta['ExposureTime'] = float(val.num) / float(val.den)
            else:
                meta['ExposureTime'] = float(val)
            break
            
    # 2. Extract shared static tags
    print(tags)
    if 'Image Make' in tags:
        meta['Make'] = str(tags['Image Make']).strip()
        
    if 'Image Model' in tags:
        meta['Model'] = str(tags['Image Model']).strip()
        
    if 'EXIF ISOSpeedRatings' in tags:
        meta['ISO'] = int(tags['EXIF ISOSpeedRatings'].values[0])
        
    if 'EXIF FNumber' in tags:
        val = tags['EXIF FNumber'].values[0]
        if hasattr(val, 'num') and hasattr(val, 'den'):
            meta['FNumber'] = (val.num, val.den)
            
    if 'EXIF FocalLengthIn35mmFilm' in tags:
        val = tags['EXIF FocalLengthIn35mmFilm'].values[0]
        if hasattr(val, 'num') and hasattr(val, 'den'):
            meta['FocalLengthIn35mmFilm'] = (val.num, val.den)
            
    if 'EXIF LensModel' in tags:
        meta['LensModel'] = str(tags['EXIF LensModel']).strip()
        
    return meta


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


def save_tiff_16bit(path, image_float01, icc_profile=None, shared_metadata=None):
    """image_float01: HxWx3 float in [0,1] -> written as 16-bit TIFF."""
    out = np.clip(image_float01, 0.0, 1.0)
    out_u16 = (out * 65535.0 + 0.5).astype(np.uint16)
    
    extratags = []
    if icc_profile is not None:
        extratags.append((ICC_TAG_CODE, "B", len(icc_profile), icc_profile, True))
        
    if shared_metadata:
        # Structure: (tag_code, type_code, count, value, writeonce)
        if 'Make' in shared_metadata:
            extratags.append((271, 's', 0, shared_metadata['Make'], True))
        if 'Model' in shared_metadata:
            extratags.append((272, 's', 0, shared_metadata['Model'], True))
        if 'LensModel' in shared_metadata:
            extratags.append((42036, 's', 0, shared_metadata['LensModel'], True))
        if 'ISO' in shared_metadata:
            extratags.append((34855, 'H', 1, shared_metadata['ISO'], True))
        if 'FNumber' in shared_metadata:
            extratags.append((33437, '2I', 1, shared_metadata['FNumber'], True))
        if 'FocalLength' in shared_metadata:
            extratags.append((37386, '2I', 1, shared_metadata['FocalLength'], True))

    tifffile.imwrite(path, out_u16, photometric="rgb", extratags=extratags)


def save_preview(path, image_float01):
    out = np.clip(image_float01, 0.0, 1.0)
    out_u8 = (out * 255.0 + 0.5).astype(np.uint8)
    cv2.imwrite(path, cv2.cvtColor(out_u8, cv2.COLOR_RGB2BGR))


# --------------------------------------------------------------------------
# Alignment
# --------------------------------------------------------------------------

def to_gray_f32(img_u16):
    f = img_u16.astype(np.float32)
    gray = 0.2126 * f[..., 0] + 0.7152 * f[..., 1] + 0.0722 * f[..., 2]
    gray /= 65535.0
    return gray


def build_pyramid(img, levels):
    pyr = [img]
    for _ in range(levels - 1):
        pyr.append(cv2.pyrDown(pyr[-1]))
    return pyr[::-1]


def ecc_align_pair(template_gray, input_gray, motion_type, levels=4,
                   iterations=300, eps=1e-7):
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
            break
        if lvl < levels - 1:
            warp[:, 2] *= 2.0
    return warp


def _to_3x3(m2x3):
    m = np.eye(3, dtype=np.float64)
    m[:2, :] = m2x3
    return m


def _to_2x3(m3x3):
    return m3x3[:2, :].astype(np.float32)


def align_stack(images, ref_index, motion_type, levels=4):
    n = len(images)
    grays = [to_gray_f32(im) for im in images]
    cumulative = [None] * n
    cumulative[ref_index] = np.eye(3, dtype=np.float64)

    for k in range(ref_index - 1, -1, -1):
        print(f"  aligning frame {k} -> {k + 1}")
        pairwise = ecc_align_pair(grays[k + 1], grays[k], motion_type, levels)
        cumulative[k] = _to_3x3(pairwise) @ cumulative[k + 1]

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
    imgs_255 = [im.astype(np.float32) * (255.0 / 65535.0) for im in images_u16]
    merger = cv2.createMergeMertens()
    fused = merger.process(imgs_255)
    return np.clip(fused, 0.0, 1.0)


def merge_debevec(images_u16, exposure_times_s, tonemap_method="reinhard", gamma=1.0):
    imgs_8u = [
        np.clip(im.astype(np.float32) / 257.0, 0, 255).astype(np.uint8)
        for im in images_u16
    ]
    times = np.array(exposure_times_s, dtype=np.float32)
    calibrate = cv2.createCalibrateDebevec()
    response = calibrate.process(imgs_8u, times)
    merger = cv2.createMergeDebevec()
    hdr = merger.process(imgs_8u, times, response)

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
                   default="translation", help="Alignment motion model")
    p.add_argument("--ref-index", type=int, default=None,
                   help="Index of the reference frame. Default: the middle frame.")
    p.add_argument("--pyramid-levels", type=int, default=4,
                   help="Number of pyramid levels for ECC alignment (default: 4)")
    p.add_argument("--max-dimension", type=int, default=None,
                   help="Optionally downscale so the longer side is at most this many pixels")
    p.add_argument("--exposures", type=float, nargs="+", default=None,
                   help="Optional manual exposure times in seconds. If omitted, "
                        "they are read automatically from the TIFF headers.")
    p.add_argument("--auto-exposures", action="store_true",
                   help="Auto-generate a theoretical exposure sequence using base/ev-step.")
    p.add_argument("--base-exposure", type=float, default=0.5,
                   help="Longest exposure in seconds, for --auto-exposures")
    p.add_argument("--ev-step", type=float, default=1.7,
                   help="EV step between brackets, for --auto-exposures")
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
        print("\nEmbedded ICC profile found in first file; carrying over to output.")

    # Extract metadata automatically
    print("\nReading EXIF metadata from files...")
    all_metadata = []
    for f in files:
        all_metadata.append(get_file_metadata(f))
        
    # Determine which metadata tags are shared across ALL images
    shared_metadata = dict(all_metadata[0])
    for meta in all_metadata[1:]:
        keys_to_remove = []
        for k in shared_metadata:
            if k not in meta or meta[k] != shared_metadata[k]:
                keys_to_remove.append(k)
        for k in keys_to_remove:
            del shared_metadata[k]
            
    # Exposure time varies intentionally, so we shouldn't burn one into the HDR metadata
    if 'ExposureTime' in shared_metadata:
        del shared_metadata['ExposureTime']
        
    print(f"Shared static metadata detected: {list(shared_metadata.keys())}")

    print("\nLoading images...")
    images = load_images(files, max_dimension=args.max_dimension)
    n = len(images)

    ref_index = args.ref_index if args.ref_index is not None else n // 2
    if not (0 <= ref_index < n):
        sys.exit(f"--ref-index {ref_index} out of range for {n} images")
    print(f"\nReference frame: index {ref_index} ({files[ref_index]})")

    # Resolve exposure times for Debevec
    exposures = None
    if args.method == "debevec":
        if args.auto_exposures:
            exposures = auto_exposure_sequence(n, args.base_exposure, args.ev_step)
            print("  using auto-generated exposure sequence (seconds):",
                  [f"{e:.5f}" for e in exposures])
        elif args.exposures is not None:
            if len(args.exposures) != n:
                sys.exit(f"--exposures has {len(args.exposures)} values but there are {n} images.")
            exposures = args.exposures
        else:
            # Fall back to EXIF reading
            exposures = []
            for i, meta in enumerate(all_metadata):
                if 'ExposureTime' in meta:
                    exposures.append(meta['ExposureTime'])
                else:
                    sys.exit(f"Could not read exposure time from {files[i]}. "
                             "Please provide --exposures manually.")
            print("  using exposure times read from EXIF (seconds):", 
                  [f"{e:.5f}" for e in exposures])

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
    save_tiff_16bit(args.output, fused, icc_profile=icc_profile, shared_metadata=shared_metadata)

    if args.preview:
        preview_path = os.path.splitext(args.output)[0] + "_preview.jpg"
        save_preview(preview_path, fused)
        print(f"Saved preview {preview_path}")

    print("Done.")


if __name__ == "__main__":
    main()
