function writeNIFTI(name, data, info)


OptDelOut.Prefix = name;
command = ['rm -f ' name '.nii.gz']; [status,cmdout] = system(command,'-echo');
[err, ErrMessage, InfoDelOut] = WriteBrik(data, info, OptDelOut);

command = ['3dAFNItoNIFTI -overwrite -prefix ' name '.nii.gz ' name '+orig'];
[status,cmdout] = system(command,'-echo');
command = ['rm ' name '+orig*'];
[status,cmdout] = system(command,'-echo');