   j = 32;
   newbiomass = 10;
   mergemodeltest40001test1=importExcelModel('/mnt/NFS/fengch/new/models/merge.xlsx')
   load ('mergemodeltest40001test1.mat')
   s = full(mergemodeltest40001test1.S);
   mergemodeltest40001test1.lb([1663 1696])=[0 -1000]
   mergemodeltest40001test1.ub([1663 1696])=[0 0]
   %mergemodeltest40001test1.c([2 3])=[1 0]
   % disp(mergemodeltest40.c(2458))
   one_point_pos = [2963:3023]';
   three_point_pos = [3024:5965]';
   length_z = length(one_point_pos)+length(three_point_pos);

   new_s = zeros(size(s,1)+2*length_z,size(s,2)+length_z);
   new_s(1:size(s,1),1:size(s,2)) = s;
   for i = 1:length(one_point_pos)
      new_s(size(s,1)+i,one_point_pos(i)) = 1;
      new_s(size(s,1)+i,size(s,2)+i) = 1000;
      new_s(size(s,1)+length_z+i,one_point_pos(i)) = 1;
      new_s(size(s,1)+length_z+i,size(s,2)+i) = -1000;
   end

   for i = 1:length(three_point_pos)
      new_s(size(s,1)+length(one_point_pos)+i,three_point_pos(i)) = 1;
      new_s(size(s,1)+length(one_point_pos)+i,size(s,2)+length(one_point_pos)+i) = 1000;
      new_s(size(s,1)+length_z+length(one_point_pos)+i,three_point_pos(i)) = 1;
      new_s(size(s,1)+length_z+length(one_point_pos)+i,size(s,2)+length(one_point_pos)+i) = -1000;
   end

   c = zeros(size(new_s,2),1);
   c(size(s,2)+1:size(s,2)+length(one_point_pos)) = 1;
   c(size(s,2)+length(one_point_pos)+1:size(s,2)+length_z) = 3;
   %c(find(mergemodeltest40.c)) = 1;

   blc = zeros(size(new_s,1),1);
   buc = zeros(size(new_s,1),1);
   buc(size(s,1)+1:size(s,1)+length_z) = 2000;
   buc(size(s,1)+length_z+1:end) = 0;
   blc(size(s,1)+length_z+1:end) = -2000 + 0.00001;

   blx = zeros(size(new_s,2),1);
   bux = zeros(size(new_s,2),1);
   bux(1:size(s,2)) = 1000;
   bux(size(s,2)+1:end) = 1;
   rev_pos = find(mergemodeltest40001test1.lb < 0);
   blx(rev_pos) = -1000;
   %{
   blx(one_point_pos) = -1000;
   blx(three_point_pos) = -1000;
   %}
   blx(find(mergemodeltest40001test1.c)) = 0.1*newbiomass;

   int_sub = [size(s,2)+1:size(s,2)+length_z];
   prob.a = new_s;
   prob.c = c;
   prob.blc = blc;
   prob.buc = buc;
   prob.blx = blx;
   prob.bux = bux;

   prob.ints.sub = int_sub;
   [r,res] = mosekopt('minimize',prob);
      result = res.sol.int.pobjval;
   resultss = find(res.sol.int.xx(size(s,2)+1:end));
   currentFile=sprintf('ggtrans%d.mat',j)
   save(currentFile,'resultss')

%%params = readtable('/mnt/NFS/fengch/new/transport_matlab_new.csv')
%%for j = 29:39
%%newex = params.Var1(j);
%%newbiomass = params.Var2(j);
%%disp(newex)