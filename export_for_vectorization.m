function []= export_for_vectorization(filename,narrowBand,theta1_output,theta2_output,tauBlurred)

narrowBand(:,end)=[];%remove padded column to the right
narrowBand(end,:)=[];%remove bottom padded row

f = fopen(filename,'w');
fprintf(f, 'narrowBand = [');
m=size(narrowBand,1);
n =size(narrowBand,2);
for i=1:m
 for j=1:n
  fprintf(f,"%d",narrowBand(i,j));
  if (j~=n)
   fprintf(f," ");
  end
 end
 fprintf(f,";\n");
end
fprintf(f,'];\n');

fprintf(f, 'theta1_output = [');
for i=1:m
 for j=1:n
  fprintf(f,"%f",theta1_output(i,j));
  if (j~=n)
   fprintf(f," ");
  end
 end
 fprintf(f,";\n");
end
fprintf(f,'];\n');

fprintf(f, 'theta2_output = [');
for i=1:m
 for j=1:n
  fprintf(f,"%f",theta2_output(i,j));
  if (j~=n)
   fprintf(f," ");
  end
 end
 fprintf(f,";\n");
end
fprintf(f,'];\n');

fprintf(f, 'tauBlurredRe = [');
for i=1:m
 for j=1:n
  fprintf(f,"%f",real(tauBlurred(i,j)));
  if (j~=n)
   fprintf(f," ");
  end
 end
 fprintf(f,";\n");
end
fprintf(f,'];\n');

fprintf(f, 'tauBlurredIm = [');
for i=1:m
 for j=1:n
  fprintf(f,"%f",imag(tauBlurred(i,j)));
  if (j~=n)
   fprintf(f," ");
  end
 end
 fprintf(f,";\n");
end
fprintf(f,'];\n');


fclose(f);

end