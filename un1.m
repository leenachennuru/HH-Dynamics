Vplotfinal = [];
for  m = 1:b
   for k = 1:t
Vplotneeded{m,k} = Vsave{m,k}(7000:10000);
[vplotx,vploty] = size(Vplotneeded{m,k});
Vfinplot{m,k} = [Vplotneeded{m,k};Iext(m)*ones(vplotx,vploty)];
    end
end

Vplotfinal = zeros(2,603201);

for m = 1:b
    Vplotfinal(:,((m-1)*1000)+1 : ((m*1000)+1)) = Vfinplot{m,2};
end

for m = 1:b
Vplotfinal = [Vplotfinal:Vfinplot{m,2}];
end


figure,plot(Vplotfinal(2,:),Vplotfinal(1,:))
Vsave cell with b rows .. and m columns ...  