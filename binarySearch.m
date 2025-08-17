function idx = binarySearch(A, num)
   l = 1;
   r = length(A);
   idx = 1;
   while l < r
      idx = 1 + floor((l + r - 1) / 2);
      if A(idx) > num
        r = idx - 1; 
      elseif A(idx) <= num
        l = idx;
      end
   end
   if l == r
      idx = r; 
   end
   if A(idx) > num
     idx = -1;
   end
end