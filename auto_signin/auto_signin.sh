host=$1
export LC_CTYPE="en_US.UTF-8"
expect -c "
spawn ssh fangzhonghai@$host
set timeout 30
expect  \"Password:\"
set password \"BGI_FZH_2022\"
send \"\$password\r\"
expect \"Verification code:\"
set key [ exec python3 /home/fangzhonghai/get_totp_code.py ]
send \"\$key\r\"
interact
"
