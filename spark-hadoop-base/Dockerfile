FROM fdiblen/archsci:latest
MAINTAINER fdiblen


USER root
#USER archsci
ENV DISPLAY :0


# INSTALL EXTRA PACKAGES
# ===================================
RUN pacman -Syyuu --noconfirm --needed
RUN pacman -S --noconfirm --needed \
    jdk8-openjdk \
    ipython2 ipython r rsync
RUN archlinux-java set java-8-openjdk/jre

RUN mkdir -p /opt/soft

ADD scripts/install_spark.sh /tmp/
RUN chmod +x /tmp/install_spark.sh && \
    bash /tmp/install_spark.sh && \
    rm -f /tmp/install_spark.sh

ADD scripts/install_hadoop.sh /tmp/
RUN chmod +x /tmp/install_hadoop.sh && \
    bash /tmp/install_hadoop.sh && \
    rm -f /tmp/install_hadoop.sh

COPY conf/spark-hadoop-env.sh  /etc/profile.d/

# SETUP SSH FOR HADOOP USER
# ===================================
#RUN sudo -H -u root bash -c 'cat /etc/hadoop/.ssh/id_rsa.pub >> /home/archsci/.ssh/authorized_keys' && \
#    sudo -H -u root bash -c 'cat /etc/hadoop/.ssh/id_rsa.pub >> /root/.ssh/authorized_key'

RUN sudo pacman --noconfirm -Scc && \
    sudo rm -rf /home/archsci/temp && \
    sudo rm -rf /home/archsci/.cache && \
    sudo rm -rf /var/cache/pacman/pkg/* \
    rm -rf /tmp \
    sudo rm -rf yaourt-tmp-*


ENV JAVA_HOME=/usr/lib/jvm/default
ENV SPARK_VERSION=2.2.1
ENV HADOOP_VERSION=2.7.5

ENV SPARK_HOME=/opt/soft/spark
ENV HADOOP_HOME=/opt/soft/hadoop

ENV HADOOP_NATIVE=/opt/soft/hadoop/lib/native/
ENV JAVA_LIBRARY_PATH=$JAVA_LIBRARY_PATH:$HADOOP_NATIVE
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HADOOP_NATIVE:/opt/soft/casacore/lib
ENV SPARK_YARN_USER_ENV="JAVA_LIBRARY_PATH=$JAVA_LIBRARY_PATH,LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

ENV HADOOP_COMMON_LIB_NATIVE_DIR=/opt/soft/hadoop/lib/native
ENV HADOOP_OPTS="$HADOOP_OPTS -Djava.library.path=$HADOOP_HOME/lib/native"

ENV SPARK_MASTER_PORT 7077
ENV SPARK_MASTER_WEBUI_PORT 8080
ENV SPARK_MASTER_LOGFILE=/opt/soft/spark/logs/spark_master.log


EXPOSE 22
#CMD ["/bin/zsh"]
CMD ["sudo","/usr/bin/supervisord", "-c", "/etc/supervisor/conf.d/supervisord.conf"]
