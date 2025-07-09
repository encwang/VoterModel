function y = NchooseK(i,j)

if i < j
    y = 0;
else
    y = nchoosek(i,j);
end