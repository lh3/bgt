EXE=./bgt

if [ ! -x $EXE ]; then
	echo "ERROR: failed to find '$EXE' executable."
	exit 1
fi

if [ ! -f 1kg11-1M.raw.bcf ] || [ ! -f 1kg11-1M.raw.samples.gz ] || [ ! -f anno11-1M.fmf.gz ]; then
	echo "MESSAGE: downloading example data..."
	wget -qO- http://bit.ly/BGTdemo | tar xf -
fi

MD5=md5sum
which md5sum
if [ $? -ne 0 ]; then
	MD5="md5 -r"
fi

echo "MESSAGE: importing..."
$EXE import 1kg11-1M.bgt 1kg11-1M.raw.bcf
gzip -dc 1kg11-1M.raw.samples.gz > 1kg11-1M.bgt.spl

echo "MESSAGE: computing checksum of various queries..."
$MD5 1kg11-1M.bgt.bcf | awk '{print $1}'
$EXE view -C 1kg11-1M.bgt | $MD5
$EXE view -s,HG00171,HG00173 -f'AC>0' -r 11:100000-200000 1kg11-1M.bgt | $MD5
$EXE view -s'population=="CEU"' -s'population=="YRI"' -f'AC1/AN1>=0.1&&AC2==0' -G 1kg11-1M.bgt | $MD5
$EXE view -d anno11-1M.fmf.gz -a'impact=="HIGH"' -CG 1kg11-1M.bgt | $MD5

echo -e "\nCorrect checksum should be:"
echo 03cc454d4495500061a82ebe341b89ef
echo 17a6a271f0057b6dbfe9b44c6317c5c7
echo 76633b2f9efe8d5b8b39868bb24a51f4
echo 722ae5f5671c4e024842c59f80a11d16
echo 7709cceaec9a1f084e3f509a72a7a615
