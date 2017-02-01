function [file_names,classes] = read_gw_txt(filename)

fp = fopen(filename) ;
S = textscan(fp,'%s %s') ;
fclose(fp) ;

file_names = S{2} ;
classes = S{1} ;

end