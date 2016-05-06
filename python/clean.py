#!/usr/bin/env python3

import sys
import re

with open(sys.argv[1], mode='r', newline='\n') as fin,\
     open(sys.argv[2], mode='w', newline='\n') as fout:
        for line in fin:
            fout.write(re.sub('.*\r', '', line))
