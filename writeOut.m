 clear all; 
 clc;

load result;

for i = 1 : 9
    filename=['MTOSOO_P' num2str(i) '.txt'];
    fid=fopen(filename,'w');
    for r = 20 : 20 : 2001
        fprintf(fid,'%d,',r*100,MFEA_TLS(i).EvBestFitness(1:59,r));
        fprintf(fid,'%d',MFEA_TLS(i).EvBestFitness(60,r));
        fprintf(fid,'\r\n');
    end
    fclose(fid);
end

