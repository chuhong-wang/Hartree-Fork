for i in {0..9}
do
cp 1.vesta $i.vesta
gsed -i "s/REPLACE/$i/g" $i.vesta
done
