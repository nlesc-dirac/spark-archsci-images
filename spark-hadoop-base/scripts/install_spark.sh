#!/bin/sh

echo ""
echo "      SPARK INSTALLATION"
echo ""

cd $(mktemp -d)
pwd

version=2.2.1

wget -c --quiet http://archive.apache.org/dist/spark/spark-$version/spark-$version-bin-without-hadoop.tgz
tar -xzf spark-$version*
rm -rf spark*.tgz
mv spark-* /opt/soft/
ln -s /opt/soft/spark-* /opt/soft/spark




# sudo -H -u archsci bash -c 'yaourt -S --needed --noconfirm apache-spark'
# ln -s /opt/apache-spark/ /opt/soft/spark
