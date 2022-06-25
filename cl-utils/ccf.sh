http://www.cosmos.esa.int/web/xmm-newton/current-calibration-files
rsync -v -a --delete --delete-after --force --include='*.CCF' --exclude='*/' xmm.esac.esa.int::XMM_VALID_CCF /home/dalp/ccf/
rsync -v -a --delete --delete-after --force --include='*.CCF' --exclude='*/' xmm.esac.esa.int::XMM_VALID_CCF /Users/silver/ccf/

latest update 2019 Dec 11
