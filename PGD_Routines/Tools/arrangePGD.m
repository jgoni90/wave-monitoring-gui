function pgd = arrangePGD(pgd1,pgd2,f1,f2,OPTION)

% function pgd = arrangePGD(pgd1,pgd2,f1,f2,OPTION)
%
% OPTION = 'add' returns pgd.RB{i} = [f1(i)*pgd1.RB{i}, f2(i)*pgd2.RB{i}] for all
% i=1,..,nsd.
% NOTE: factors f1 and f2 have nsd components.
%
% OPTION = 'merge' returns pgd.RB{i}(:,j) = f1*pgd1.RB{i}(:,j) + f2*pgd2.RB{i}(:,j)
% for all i=1,..,nsd and j=1,..,nmodes even if nmodes is not the same for pgd1 and pgd2
% NOTE: factors f1 and f2 are scalar.
%
% OPTION = 'prod' returns pgd.RB{i}(:,k) = pgd1.RB{i}(:,j) .* pgd2.RB{i}(:,p)
% for all i=1,..,nsd, j=1,...,nmodes1, p=1,...,nmodes2 and k=1,..,nmodesl*nmodes2
% NOTE1: input f1(i)=true returns sparse structure pgd.RB{i}(f2{i}{1},k)
% with size f2{i}{2} x nmodesl*nmodes2.
% NOTE2: input f1(i)=false ignores f2{i}.
%
% OPTION = 'extend' returns pgd.RB{i}(f2{i}{1},j) = pgd1.RB{i}(:,j) for all 
% i=1,..,nsd, j=1,...,nmodes1, with size f2{i}{2} x nmodesl, if and only if f1(i)=true.
% NOTE1: input pgd2 is ignored.
% NOTE2: input f1(i)=false ignores f2{i} and it returns pgd.RB{i} = pgd1.RB{i}.
%
% OPTION = 'extends' returns sparse structure pgd.RB{i} of OPTION = 'extend'.
%
% OPTION = 'join2' returns pgd.RB{f1}(:,j) = pgd1.RB{f1}(:,j) x pgd1.RB{f2}(:,j).' 
% for all j.
% NOTE1: pgd2 is empty and f1, f2 are scalars.
% NOTE2: pgd.RB(f2) will be deleted.

switch OPTION
    case 'merge'
        [maxterms,posmax] = max([pgd1.counters.Uterm,pgd2.counters.Uterm]);

        if posmax == 1
            pgd = pgd1;
            for j = 1:numel(pgd.RB)
                pos1 = 1:pgd2.counters.Uterm;
                pos2 = pgd2.counters.Uterm+1:maxterms;
                pgd.RB{j}(:,pos1) = f1*pgd1.RB{j}(:,pos1) + f2*pgd2.RB{j}(:,pos1);
                pgd.RB{j}(:,pos2) = f1*pgd1.RB{j}(:,pos2) + f2*0;
            end
        else
            pgd = pgd2;
            for j = 1:numel(pgd.RB)
                pos1 = 1:pgd1.counters.Uterm;
                pos2 = pgd1.counters.Uterm+1:maxterms;
                pgd.RB{j}(:,pos1) = f1*pgd1.RB{j}(:,pos1) + f2*pgd2.RB{j}(:,pos1);
                pgd.RB{j}(:,pos2) = f1*0 + f2*pgd2.RB{j}(:,pos1);
            end
        end
    
    case 'add'
        pgd = pgd1;
        for j = 1:numel(pgd.RB)
            pgd.RB{j} = [f1(j)*pgd1.RB{j}, f2(j)*pgd2.RB{j}];
        end
        pgd.counters.Uterm = size(pgd.RB{1},2);
        
    case 'prod'
        pgd.RB = cell(numel(pgd1.RB),1);
        
        for j = 1:numel(pgd1.RB)     
            jnodes1 = size(pgd1.RB{j},1);
            jterms1 = size(pgd1.RB{j},2);
            jterms2 = size(pgd2.RB{j},2);
            
            if f1(j) %Sparse structure
                
               termsv = 1:jterms1*jterms2;
               I = f2{j}{1}(:,ones(size(termsv)));
               J = termsv(ones(size(f2{j}{1})),:);
               V = zeros(size(I));
               index = 1;
                for i = 1:jterms1
                    for k = 1:jterms2
                        V(:,index) = pgd1.RB{j}(:,i) .* pgd2.RB{j}(:,k);
                        index = index + 1;
                    end
                end
                pgd.RB{j} = sparse(I(:),J(:),V(:),f2{j}{2},jterms1*jterms2);
                
            else
                
                pgd.RB{j} = zeros(jnodes1,jterms1*jterms2);
                index = 1;
                for i = 1:jterms1
                    for k = 1:jterms2
                        pgd.RB{j}(:,index) = pgd1.RB{j}(:,i) .* pgd2.RB{j}(:,k);
                        index = index + 1;
                    end
                end
                
            end
        end
        
    case 'extend'
        pgd = pgd1;
        
        for j = 1:numel(pgd.RB)
            if f1(j) %Do extension
                jterms1 = size(pgd1.RB{j},2);
                pgd.RB{j} = zeros(f2{j}{2},jterms1);
                pgd.RB{j}(f2{j}{1},:) = pgd1.RB{j};
            end
        end
        
    case 'extends'
        pgd = pgd1;

        for j = 1:numel(pgd.RB)
            if f1(j) %Do extension
                jterms1 = size(pgd1.RB{j},2);
                termsv = 1:jterms1;
                I = f2{j}{1}(:,ones(size(termsv)));
                J = termsv(ones(size(f2{j}{1})),:);
                pgd.RB{j} = sparse(I(:),J(:),pgd1.RB{j}(:),f2{j}{2},jterms1);
            end
        end
        
    case 'join2'
        pgd = pgd1;
        
        [n1,nt] = size(pgd.RB{f1});
        n2 = size(pgd.RB{f2},1);
        pgd.RB{f1} = zeros(n1*n2,nt);
        for j = 1:nt
            mode2 = pgd1.RB{f1}(:,j) * pgd1.RB{f2}(:,j).';
            pgd.RB{f1}(:,j) = mode2(:);
        end
        pgd.RB(f2) = [];
end






