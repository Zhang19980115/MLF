function TimeSeries=GMRecord_Shaper(Groundmotion_record,Record_source)

GM_info = importdata('8GMinfo.txt');
isMatch = ismember(Groundmotion_record, GM_info.textdata);
if isMatch==0
    disp('Cant find matched ground motion');
        keyboard
else 
   index   = find(ismember(GM_info.textdata, Groundmotion_record));
end 

dt=GM_info.data(index,1);
Numstep=GM_info.data(index,2);
t=linspace(0, (Numstep-1)*dt, Numstep)';

ug = importdata([Record_source Groundmotion_record '.txt']);

ug=ug/max(abs(ug)); % normalization
TimeSeries=[t,ug];
end 