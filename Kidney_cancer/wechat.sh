xhost +
docker run \
    -i \
    --rm \
    --ipc=host \
    --privileged \
    -e DISPLAY=$DISPLAY \
    -e XMODIFIERS=@im=fcitx \
    -e QT_IM_MODULE=fcitx \
    -e GTK_IM_MODULE=fcitx \
    -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
    -v $HOME/.config/weixin:/root/.config/weixin \
    wechat:0.0.1
xhost -