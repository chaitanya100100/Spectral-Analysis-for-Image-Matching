% this script run `main` for entire dataset.

li = dir('../dataset') ;
li = li(3:size(li)) ;

for i = 1:size(li)
    if li(i).isdir
        disp(i);
        disp(li(i).name);
        main(li(i).name);
    end
end