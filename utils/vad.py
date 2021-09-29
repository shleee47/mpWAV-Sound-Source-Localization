import numpy as np
import pdb
import librosa
import os 
import soundfile
import random
import pickle
import torch
import torchaudio
import math

wav_path = "../pickle/1_enhanced/"
out_path = "../pickle/2_vad/"
os.makedirs(out_path,exist_ok=True)

wav_list = [f for f in os.listdir(wav_path)]
wav_list = sorted(wav_list)
frame_size = 512
hop_size = 128
RISING_TERM = 30 # 15*4
LEAST_TERM = 80 #25*4
power_list = []
save_wav_count=0
eps = 1e-5

def stft(_wav_path):
    window=torch.hann_window(window_length=512, periodic=True, dtype=None, layout=torch.strided, device=None, requires_grad=False)
    data_wav,_ = librosa.load(_wav_path,sr=16000,mono=False)
    #pdb.set_trace()
    data_wav = torch.from_numpy(data_wav)
    spec_noi1 = torchaudio.functional.spectrogram(waveform=data_wav, pad=0, window=window, n_fft=512, hop_length=128, win_length=512, power=None, normalized=False)
    input_wav_real1 =spec_noi1[:,:,:,0]
    input_wav_imag1 = spec_noi1[:,:,:,1]
    phase = torch.atan(input_wav_imag1/(input_wav_real1+1e-8))
    input_wav_magnitude = torch.sqrt(input_wav_real1**2 + input_wav_imag1**2)
    return input_wav_magnitude, phase

#def spec2audio_tensor(power_mag,phase,window,length,nfft):
#    window = window
#    length = length
#    mag = power_mag  #[1 F T]
#    phase = phase #[1 F T]
#    sqrt_mag = torch.sqrt(mag)
#    cos_phase = torch.cos(phase)
#    sin_phase = torch.sin(phase)
#    real = sqrt_mag * cos_phase
#    imagine = sqrt_mag * sin_phase
#    real = real.unsqueeze(3)
#    imagine = imagine.unsqueeze(3)
#    complex_ri = torch.cat((real,imagine),3)
#    audio = torch.istft(input = complex_ri, n_fft=int(nfft), hop_length=int(0.25*nfft), win_length=int(nfft), window=window, center=True, normalized=False, onesided=True, length=length)
#    return audio

#pdb.set_trace()
for pkl in wav_list:
    time_list=[]
    left_or_right_list=[]
    
    with open(os.path.join(wav_path,pkl),'rb') as f:
        data = pickle.load(f)
    wav = data["output_path"]
    #wav = data["audio_path"]

    power_list=[]
    active_cnt = 0
    tmp_dummy_frame = 0
    dummy_frame = 0
    inactive_cnt = 0
    state= 0
    i=0
    num=0 
    #pdb.set_trace()
    (total_audio,fs) = soundfile.read(wav)    
    if len(total_audio.shape) == 1:
        pdb.set_trace()
    if fs != 16000:
        pdb.set_trace()
        tmp0 = librosa.resample(total_audio[:,0],fs,16000)
        tmp1 = librosa.resample(total_audio[:,1],fs,16000)
        total_audio = np.stack((tmp1,tmp2),axis=1)
        fs = 16000
    else:
        tmp0 = total_audio[:,0]
        tmp1 = total_audio[:,1]
    
    
    frame_idx_list = range(0,len(total_audio)-hop_size+1,hop_size)
    input_wav_mag,phase = stft(wav)
        
    mean_power = abs(input_wav_mag[:,20:,:]).mean()
    thre = mean_power / 10
    
 
    for frame_idx in frame_idx_list:
        num+=1
        if abs(input_wav_mag[:,20:,frame_idx//hop_size]).mean() > thre:
            if state == 0:
                active_cnt = 1
                tmp_dummy_frame = 1
                rising_idx = frame_idx
                state =1

            elif state == 1:
                active_cnt+=1
                tmp_dummy_frame+=1
                if active_cnt == RISING_TERM:
                    state=2

            elif state == 2:
                active_cnt+=1

            elif state == 3:
                inactive_cnt=0
                active_cnt+=1
                state = 2
            
            elif state == 4:
                active_cnt =1
                tmp_dummy_frame = 1
                rising_idx = frame_idx
                state = 1

        else:

            if state == 0:
                dummy_frame+=1
                state = 0

            elif state == 1:
                active_cnt = 0
                dummy_frame+=tmp_dummy_frame
                tmp_dummy_frame = 0
                state=0
            
            elif state == 2:
                inactive_cnt =1
                active_cnt+=1
                state = 3

            elif state == 3:
                inactive_cnt+=1
                active_cnt+=1
                if inactive_cnt == LEAST_TERM:
                    state = 4

            elif state == 4:
                dummy_frame = 1
                state = 0 

        # save VAD chunk here in wav
        if state == 4 or (num == len(frame_idx_list) and active_cnt > RISING_TERM):
            falling_idx = frame_idx
            if rising_idx-hop_size < 0:
                rising_idx = 128
            rising_idx = (rising_idx-hop_size)
            if state == 4:
                falling_idx = (falling_idx-(LEAST_TERM-2)*hop_size)
            else:
                falling_idx = (falling_idx-(inactive_cnt-2)*hop_size)
            tmp0_power = np.sum(np.abs(tmp0[rising_idx:falling_idx]))
            tmp1_power = np.sum(np.abs(tmp1[rising_idx:falling_idx]))
            if tmp0_power > tmp1_power:
                left_or_right_list.append(0)
            else:
                left_or_right_list.append(1)

            rising_idx = rising_idx/fs
            falling_idx = falling_idx/fs
            time_list.append([rising_idx,falling_idx])
            save_wav_count +=1
            #save chunk for another channel 
            i+=1
            state = 4
            active_cnt = 0
            inactive_cnt = 0
            tmp_dummy_frame = 0
            dummy_frame = 0

    #pdb.set_trace()
    data["time"] = time_list
    data["LR"] = left_or_right_list
    wav_total_mean_power = np.mean(np.abs((tmp0+tmp1)/2))
    if time_list == [] or wav_total_mean_power < eps:
        pdb.set_trace()
        wav = data["input_path"]
        (total_audio,fs) = soundfile.read(wav)    
        
        if fs != 16000:
            tmp0 = librosa.resample(total_audio[:,0],fs,16000)
            tmp1 = librosa.resample(total_audio[:,1],fs,16000)
            total_audio = np.stack((tmp1,tmp2),axis=1)
            fs = 16000
        else:
            tmp0 = total_audio[:,0]
            tmp1 = total_audio[:,1]

        tmp0_power = np.sum(np.abs(tmp0))
        tmp1_power = np.sum(np.abs(tmp1))
        if tmp0_power > tmp1_power:
            left_or_right_list.append(0)
        else:
            left_or_right_list.append(1)
    
        data["LR"] = left_or_right_list
    
    with open(os.path.join(out_path,pkl),"wb") as fw:
        pickle.dump(data,fw)
    print("pickle dumped!!: {}".format(pkl))