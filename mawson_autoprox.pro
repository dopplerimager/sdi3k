flis = findfile("D:\users\SDI3000\Data\Mawson\Incoming_Zips\*.zip")

for j=0,n_elements(flis) - 1 do begin
;    spawn, "c:\progra~1\7-zip\7z e -y " + flis(j)
;    spawn, "move *_stripped_* D:\users\SDI3000\Data\Mawson"

    fprse = mc_fileparse(flis(j))
    spawn, "rename " + flis(j) + " "  + fprse.name_only + '.unzipped'
endfor
end