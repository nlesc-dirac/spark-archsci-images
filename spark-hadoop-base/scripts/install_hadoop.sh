#!/bin/sh

echo ""
echo "      HADOOP INSTALLATION"
echo ""

cd $(mktemp -d)
pwd


version=2.7.6

# Install Hadoop

#wget -c --quiet http://ftp.nluug.nl/internet/apache/hadoop/common/hadoop-$version/hadoop-$version.tar.gz
wget -c --quiet http://www-eu.apache.org/dist/hadoop/common/hadoop-$version/hadoop-$version.tar.gz
tar -xzf hadoop-*.tar.gz
rm -rf hadoop*.tar.gz
mv hadoop-* /opt/soft/
ln -s /opt/soft/hadoop-* /opt/soft/hadoop

echo "export JAVA_HOME=/usr/lib/jvm/default" >> /opt/soft/hadoop/etc/hadoop/hadoop-env.sh




# sudo -H -u archsci bash -c 'yaourt -S --needed --noconfirm hadoop'
# ln -s /usr/lib/hadoop-* /opt/soft/hadoop
