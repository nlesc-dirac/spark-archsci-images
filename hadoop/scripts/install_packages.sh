#!/bin/bash

#cp /etc/pacman.d/mirrorlist /etc/pacman.d/mirrorlist.backup
#rankmirrors -n 6 /etc/pacman.d/mirrorlist.backup > /etc/pacman.d/mirrorlist

# DEPENDENCIES
#pacman -Syyuu --needed --noconfirm --force

pacman -S --needed --noconfirm \
    rsync

#sudo -H -u archsci bash -c 'yaourt -S --needed --noconfirm '

rm -rf /var/cache/pacman/pkg/* /home/archsci/temp
