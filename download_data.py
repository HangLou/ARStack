import gdown
import shutil
import os
url = "https://drive.google.com/file/d/1veYEK0gc4GNmq2mgIGWqjX0xPw_-vQMe/view?usp=share_link"
output = "data.zip"
gdown.download(url=url, output=output, fuzzy=True, quiet=False)
shutil.unpack_archive('data.zip')
os.remove('data.zip')
