# MicrobiomeHub

## shiny-server 

### 重新定义R版本
/etc/systemd/system/shiny-server.service中添加：  
Environment="R=/home/test1/miniconda3/envs/microbiome.net/bin/R"  #30site  
sudo systemctl stop shiny-server  
sudo systemctl daemon-reload   
sudo systemctl start shiny-server  

### 更改工作目录
sudo vim /etc/shiny-server/shiny-server.conf   
sudo systemctl restart shiny-server   
sudo systemctl restart shiny-server.service   


sudo systemctl stop shiny-server  
sudo systemctl daemon-reload   
sudo systemctl start shiny-server  
sudo systemctl restart shiny-server   
sudo systemctl restart shiny-server.service   
