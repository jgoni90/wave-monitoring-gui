%Script to update the file 'file' with the counter 'MSL_SAVECOUNTER'

outputfile = [file '_MSL' num2str(MSL_SAVECOUNTER)];
MSL_SAVECOUNTER = MSL_SAVECOUNTER + 1;
disp(['Saving file ' outputfile '...'])