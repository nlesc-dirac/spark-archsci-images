FROM fdiblen/spark-hadoop-base:latest
MAINTAINER fdiblen

USER root
ENV DISPLAY :0

ADD scripts/install_packages.sh /tmp/
RUN chmod +x /tmp/install_packages.sh && \
    bash /tmp/install_packages.sh && \
    rm -f /tmp/install_packages.sh

ENV HADOOP_VERSION  2.7.5
ENV HADOOP_HOME /opt/soft/hadoop
ENV HADOOP_CONF_DIR /opt/soft/hadoop/etc/hadoop/
ENV PATH $PATH:$HADOOP_HOME/bin:$HADOOP_HOME/sbin

RUN mkdir -p /opt/soft/hadoop/logs
COPY conf  $HADOOP_HOME/etc/hadoop/
RUN hdfs namenode -format
VOLUME /user

RUN mkdir -p /opt/
COPY data /opt/

COPY start_hadoop.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/*


COPY hadoop_service.sh /tmp/
RUN cat /tmp/hadoop_service.sh >> /etc/supervisor/conf.d/supervisord.conf && \
    rm -f /tmp/hadoop_service.sh


RUN sudo pacman --noconfirm -Scc && \
    sudo rm -rf /home/archsci/temp && \
    sudo rm -rf /home/archsci/.cache && \
    sudo rm -rf /var/cache/pacman/pkg/* \
    sudo rm -rf yaourt-tmp-*


EXPOSE 22 8088 9870 50030 9000 50070 50010 50020 50075 50090


CMD ["sudo","/usr/bin/supervisord", "-c", "/etc/supervisor/conf.d/supervisord.conf"]
