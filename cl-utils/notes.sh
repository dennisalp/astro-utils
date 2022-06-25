find . -type f -name "*.tar.gz" -execdir tar -xvzf {} \;
find . -name '*.gz' -exec gunzip '{}' \;

rsync -v -a -u --delete --force  /home/dalp/Dropbox/ /media/dalp/F05A1E4C5A1E1048/Dropbox/
rsync -v -a -u --delete --force  /home/dalp/ /media/dalp/F05A1E4C5A1E1048/red/

rsync -v -a --delete --delete-after --force --include='*.CCF' --exclude='*/' sasdev-xmm.esac.esa.int::XMM_VALID_CCF ~/ccf/

curl -sS -O -J "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?level=PPS&obsno=0770380401"

ffmpeg -i in.m4a out.mp3
sox --norm=-10 in.mp3 out.mp3

for f in *.jpeg; do convert -compress jpeg ./"$f" ./"${f%.jpeg}.pdf"; done

pdftk inputs*.pdf cat output output.pdf

decrypt_data.pl -d 40501004004

echo 'Message' | mail -s 'Subject' example@gmail.com
notify-send -u normal 'Topic' 'Message'

find . -name '*.reg' -maxdepth 2 | cpio -pdm ~/dat/nus/87a 

# https://dwheeler.com/sloccount/
sloccount .
