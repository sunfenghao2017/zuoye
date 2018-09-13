# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
 """
from hcacn.core import Configure, Schedule
from hcacn.steps import edgeR

#Configure.setIdentity('zywang')
Configure.enableDocker(False)
edgeR(matrixdata = "./minidata/edgeR/couts.txt", annotation = "./minidata/edgeR/test.mgrp", outputpath = None)
Schedule.run()
