"""
Copyright (c) 2014 NavPy Developers. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in
LICENSE.txt
"""

import struct
import numpy as _np

###############################################################################
###################                 READ GPS               ####################
###############################################################################
def read_gps(fid,sat):
    state = 0
    bytes_cnt = 0

    weekNo = 0
    GPS_TOW = 0
    msg_id = 0
    msgBuff = 0
    msgValid = False

    while(1):
        if(state==0):
            dataBuff = fid.read(1)
            if(dataBuff ==''):
                break
            if(dataBuff.encode('hex')=='aa'):
                Buff = dataBuff
                state = 1
                bytes_cnt += 1
    
        elif(state==1):
            dataBuff = fid.read(1)
            if(dataBuff ==''):
                break
            if(dataBuff.encode('hex')=='44'):
                Buff += dataBuff
                state = 2
                bytes_cnt += 1
            else:
                state = 0
                bytes_cnt = 0
    
        elif(state==2):
            dataBuff = fid.read(1)
            if(dataBuff ==''):
                break
            if(dataBuff.encode('hex')=='12'):
                Buff += dataBuff
                state = 3
                bytes_cnt += 1
            else:
                state = 0
                bytes_cnt = 0

        elif(state==3):
            dataBuff = fid.read(1)
            if(dataBuff ==''):
                break
            # THIS IS THE HEADER LENGTH
            hdr_len = struct.unpack('@B',dataBuff)[0]
        
            Buff += dataBuff
            state = 4
            bytes_cnt += 1

        elif(state==4):
            dataBuff = fid.read(hdr_len-bytes_cnt)
            if(dataBuff ==''):
                break
            header_data = struct.unpack('@HcBHHBBHiIHH',dataBuff)
            # UNPACK HEADER 
            # ... MESSAGE ID
            msg_id = header_data[0]
            # ... MESSAGE LENGTH
            msg_len = header_data[3]
            # ... WEEK NUMBER
            weekNo = header_data[7]
            # ... TOW
            GPS_TOW = float(header_data[8]) / 1000.0
            #print("Week Number %d, TOW = %f" %(weekNo,GPS_TOW))
        
            Buff += dataBuff
            state = 5
            bytes_cnt = 0
    
        elif(state==5):
            msgBuff = fid.read(msg_len)
            if(msgBuff == ''):
                break
            
            Buff += msgBuff
            state = 6
            bytes_cnt = 0
    
        elif(state==6):
            dataBuff = fid.read(4)
            if(dataBuff == ''):
                break
            #print dataBuff.encode('hex')
            crc_read = struct.unpack('@I',dataBuff)[0]
            crc = _CALCULATEBLOCKCRC32(Buff)
            if(crc_read != crc):
                weekNo = 0
                GPS_TOW = 0
                msg_id = 0
                msgBuff = 0
                msgValid = False
                print("TOW = %f, Bad CRC" % GPS_TOW)
            else:
                if(msg_id==43):
                    sat.set_TOW(GPS_TOW)
                    parse_range(msgBuff,sat)
                msgValid = True
                
            state = 0
            bytes_cnt = 0
            
            # NOW WE ARE DONE....
            break;
            
        else:
            break;

    return msg_id, msgValid

###############################################################################
###################                PARSE GPS               ####################
###############################################################################

def parse_range(msgBuff,sat):
    sat.set_pseudorange(_np.nan*_np.ones(32),_np.nan*_np.ones(32),range(32))
    sat.set_carrierphase(_np.nan*_np.ones(32),_np.nan*_np.ones(32),_np.nan*_np.ones(32),range(32))
    sat.set_doppler(_np.nan*_np.ones(32),range(32))
    sat.set_CNo(_np.nan*_np.ones(32),range(32))
    
    # SET IT TO THE DATA
    num_sv = struct.unpack('@I',msgBuff[0:4])[0]
    for cnt in range(num_sv):
        sv = struct.unpack('@H',msgBuff[4+cnt*44:6+cnt*44])[0]-1
        if(sv>32):
            continue
        
        #print sv
        
        sat.set_pseudorange(struct.unpack('@d',msgBuff[8+cnt*44:16+cnt*44])[0],\
                            struct.unpack('@f',msgBuff[16+cnt*44:20+cnt*44])[0],\
                            sv)
        sat.set_carrierphase(struct.unpack('@d',msgBuff[20+cnt*44:28+cnt*44])[0],\
                            struct.unpack('@f',msgBuff[28+cnt*44:32+cnt*44])[0],\
                            struct.unpack('@f',msgBuff[40+cnt*44:44+cnt*44])[0],\
                            sv)
        sat.set_doppler(struct.unpack('@f',msgBuff[32+cnt*44:36+cnt*44])[0],sv)
        sat.set_CNo(struct.unpack('@f',msgBuff[36+cnt*44:40+cnt*44])[0],sv)
        
        sat.check_dataValid(sv)

###############################################################################

###############################################################################
###################              CRC32 FUNTIONS            ####################
###############################################################################
def _CRC32VALUE(i):
    crc = i
    for j in range(8,0,-1):
        if(crc & 1):
            crc = (crc>>1) ^ 0XEDB88320
        else:
            crc = crc>>1
    return crc
    
def _CALCULATEBLOCKCRC32(entire_buffer):
    buf_len = len(entire_buffer)
    
    crc = 0
    temp1 = 0
    cnt = 0
    while(buf_len>0):
        #print crc
        temp1 = (crc >> 8) & 0X00FFFFFF
        temp2 = _CRC32VALUE( (crc ^ struct.unpack('@B',entire_buffer[cnt])[0])& 0xFF )
        
        crc = temp1 ^ temp2
        
        #print("temp1 = %d, temp2 = %d, crc = %d" %(temp1,temp2,crc))

        buf_len -= 1
        cnt += 1
    
    return crc
    
""""
import gnss

fid = open('test_files/LOG00159.BIN','rb')
satData = gnss.prn_class()

while(1):
    file_pos = fid.tell()
    
    msgID, msgValid = read_gps(fid,satData)
    
    if(msgValid):
        print("TOW = %f, msgID = %d" % (satData.get_TOW(),msgID))
    
    if( (fid.tell() - file_pos) <= 0 ):
        break

fid.close()
"""