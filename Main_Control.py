# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 17:59:32 2023
Main Python control program to interface with the Siglent SDG1032X, SDG2042X and other sililar function generators.
Calls and opens a matlab function to communicate with Siglent SDS1102X-E, SDS1104X-E

Dependency on sdg1032x.sdg1032x library https://github.com/tspspi/pysdg1032x
Dependance on matlab.engine with files specified
@author: steve
"""

import matlab.engine
import matplotlib.pyplot as plt
from sdg1032x.sdg1032x import SDG1032X
import time
import pyvisa as visa

eng = matlab.engine.start_matlab()
eng.cd(r'C:\Users\steve\Documents\School\2023Summer\Neural\KWW_Programs', nargout=0)


'Variables'
Amplitude_Jump = .1
Voltage_Start = -1
Voltage_End = 1
Frequency = 1
Port = 1
Time_Wait = 10 #seconds


'Outputs'
Ohms_C = [] 
Ohms_R = [] 
Randles_C = [] 
Randles_R = [] 
KWW_C = [] 
KWW_R = [] 
V_graph = []
KWW_N = []


'PreProcessing'
Voltage = Voltage_Start + (.5*Amplitude_Jump)
Voltage_End = Voltage_End - (.5*Amplitude_Jump) 


'Functions'
def appendAll(Ohms_CT, Ohms_RT, Randles_CT, Randles_RT, KWW_CT, KWW_RT, Voltage, KWW_NT):
    Ohms_C.append(Ohms_CT)
    Ohms_R.append(Ohms_RT)
    Randles_C.append(Randles_CT)
    Randles_R.append(Randles_RT)
    KWW_C.append(KWW_CT)
    KWW_R.append(KWW_RT)
    V_graph.append(Voltage)
    KWW_N.append(KWW_NT)

def addR2(Voltage, graph, R2):
    j = 0
    for x,y in zip(Voltage,graph):
        axs[0].annotate(float("{:.3f}".format(R2[j]), size = 6), # this is the text
                     (x,y), # these are the coordinates to position the label
                     textcoords="offset points", # how to position the text
                     xytext=(0,3), # distance from text to points (x,y)
                     ha='center') # horizontal alignment can be left, right or center
        j=j+1

def smart_query(instr, command):
    try:
        return instr.query(command)
    except visa.errors.VisaIOError:
        return None
    
def get_instr(instr_id):
    rm = visa.ResourceManager()
    if len(rm.list_resources()) == 0:
        return None
    try:
        instr = rm.open_resource(instr_id)
        return instr
    except:
        print("No instrument detected!")
        return None

'Main'
Scope = get_instr('TCPIP0::192.168.1.3::inst0::INSTR') #Oscilloscope
with SDG1032X("192.168.1.4") as sdg:
    sdg.setWaveType('SQUARE', channel=Port)
    sdg.setWaveFrequency(Frequency, channel=Port)
    sdg.setWaveOffset(Voltage, channel=Port)
    sdg.setWaveAmplitude(Amplitude_Jump,channel=Port)
    sdg.outputEnable(channel=Port, polarity="POLARITY_NORMAL", load="HZ")
    time.sleep(Time_Wait*3) #30 second wait for equalibrium
    
    while Voltage <= (Voltage_End + Amplitude_Jump):
        sdg.setWaveOffset(Voltage, channel=Port)
        TriggerLevel = (Voltage+1)*.065
        smart_query(Scope, f"C1:TRLV {TriggerLevel}")
        time.sleep(Time_Wait) # 10 second wait
        Ohms_CT, Ohms_RT, Randles_CT, Randles_RT, KWW_CT, KWW_RT, KWW_NT = eng.KWW_Record_And_Fit(Amplitude_Jump,(Voltage - (Amplitude_Jump/2)), nargout=7); #Calls matlab
        appendAll(Ohms_CT, Ohms_RT, Randles_CT, Randles_RT, KWW_CT, KWW_RT, Voltage, KWW_NT)
        
        Voltage = Voltage + Amplitude_Jump
        
    sdg.outputDisable(channel=1)


'Plotting'
fig, axs = plt.subplots(1, 2)

axs[0].plot(V_graph, Ohms_C, Label = 'Ohms Law')
axs[0].plot(V_graph, Randles_C, Label = 'Randles Fit')
axs[0].plot(V_graph, KWW_C, Label = 'KWW Fit')
axs[0].title('Effective Capacitance Vs. Voltage Offset')
axs[0].xlabel('Voltage(V)')
axs[0].ylabel('Capacitance(F)')
axs[0].legend(loc='upper right')
addR2(V_graph, KWW_C, KWW_R)

axs[1].plot(V_graph, KWW_N)
axs[1].title('N Ideality Value Vs. Voltage Offset')
axs[1].xlabel('Voltage(V)')
axs[1].ylabel('Ideality Factor')
axs[1].show()
