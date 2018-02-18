function val = mystrcmp(LIST1,LIST2,flag)

%val = mystrcmp(list1,list2,flag) Looks if all the strings of list1 are included
%in list2, no matter what order/length list1 or list2 have. An ending mark '*' 
%in any of the strings in list1 searches in list2 for the pattern 
%behind the mark. Returns boolean vector.
%
%Exception: if flag is given (optional) and activated, then mystrcmp returns false when 
%list1 and list2 have not the same length.

list1 = LIST1;
list2 = LIST2;
iscell1 = iscell(list1);
iscell2 = iscell(list2);
if iscell1, n1 = length(list1); else list1 = {list1}; n1 = 1; end
if iscell2, n2 = length(list2); else list2 = {list2}; n2 = 1; end

if nargin == 3 && flag && n1 ~= n2
    val = false;
else
    val = false(1,n1);
    for i = 1:n1
        if strcmp(list1{i}(end),'*')
            nc = length(list1{i}(1:end-1));
            for j = 1:n2, try val(i) = val(i) || strcmp(list1{i}(1:nc),list2{j}(1:nc));
                catch, val(i) = val(i) || false; end
            end
        else val(i) = any(strcmp(list1(i),list2)); 
        end
    end
end

            
            
            