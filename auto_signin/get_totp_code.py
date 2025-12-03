# -*- coding: utf-8 -*-
import pyotp
totp = pyotp.TOTP('VWBGY4IISM2ZBC2IRE6QDMB44A')
print(totp.now())
