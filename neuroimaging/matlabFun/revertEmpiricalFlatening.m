function revertEmpiricalFlatening(ax,voxProp,cont)

trans = voxProp.eccTrans.toOrig{ismember(voxProp.eccTrans.info,'empirical')};
for i = 1:length(ax.Children)
    switch ax.Children(i).Type
        case {'line' 'scatter'}
            [theta,rho] = cart2pol(ax.Children(i).XData,ax.Children(i).YData);
            [ax.Children(i).XData,ax.Children(i).YData] = pol2cart(theta,trans(rho)');
        case 'polygon'
%             [theta,rho] = cart2pol(ax.Children(i).Shape.Vertices(:,1),ax.Children(i).Shape.Vertices(:,2));
%             [X,Y] = pol2cart(theta,trans(rho));
            pgon = regions(ax.Children(i).Shape);
            for ii = 1:length(pgon)
                vert = pgon(ii).Vertices;
                [theta,rho] = cart2pol(vert(:,1),vert(:,2));
                rhot = trans(rho);
                rhot(rho<1) = min(rhot);
                [X,Y] = pol2cart(theta,rhot);
                pgon(ii).Vertices = [X Y];
                pgon(ii) = simplify(pgon(ii));
            end
            ax.Children(i).Shape = union(pgon);
        case 'image'
            im = ax.Children(i).CData;
            
            [theta,rho] = cart2pol(cont.U(cont.ind),cont.V(cont.ind));
            [U,V] = pol2cart(theta,trans(rho));
            vecUV = cont.vecUV(cont.ind);
            [X,Y] = meshgrid(ax.Children(i).XData',ax.Children(i).YData);
            
            F = scatteredInterpolant(U,V,ones(size(U)).*pi/2,'nearest','none');
            outUV = isnan(F(X,Y));
            
            F = scatteredInterpolant(U,V,vecUV,ax.Parent.UserData.F.Method,ax.Parent.UserData.F.ExtrapolationMethod);
            vecXY = F(X,Y);
            vecXY(isnan(vecXY)) = pi/2;
            
            ax.Children(i).CData = vecXY;
            
%             [theta,rho] = cart2pol(X,Y);
%             rhot = rho; rhot(:) = trans(rho);
%             [Xt,Yt] = pol2cart(theta,rhot);
%             
%             outXY = griddata(Xt,Yt,double(cont.outXY),X,Y);
%             
%             
%             cont.outXY
%             cont.X
%             cont.Y
%             
%             ind = ax.Parent.UserData.ind;
%             U = ax.Parent.UserData.U;
%             V = ax.Parent.UserData.V;
%             vecUV = ax.Parent.UserData.vecUV;
%             [theta,rho] = cart2pol(U,V);
%             [Ut,Vt] = pol2cart(theta,trans(rho));
%             [X,Y] = meshgrid(ax.Children(i).XData',ax.Children(i).YData);
%             
%             F = scatteredInterpolant(Ut(ind),Vt(ind),ones(size(Ut(ind))).*pi/2,'nearest','none');
%             outUVt = isnan(F(X,Y));
%             
%             F = scatteredInterpolant(Ut(ind),Vt(ind),vecUV(ind),ax.Parent.UserData.F.Method,ax.Parent.UserData.F.ExtrapolationMethod);
%             vecXY = F(X,Y);
%             
%             vecXY(outXY|outUV) = nan;
%             
%             
%             
%             ax.Parent.UserData.vecXY
%             ax.Children(i).CData = F(X,Y);
%             
%             
%             
%             
%             im = ax.Parent.UserData.F(Ut,Vt);
%             
%             [ax.Children(i).XData,ax.Children(i).YData] = pol2cart(theta,trans(rho)');
%             
%             ax.Parent.UserData.V
%             
%             
%             im = ax.Children(i).CData;
%             [imX,imY] = meshgrid(ax.Children(i).XData,ax.Children(i).YData');
%             
%             [theta,rho] = cart2pol(imX,imY);
%             rhot = rho;
%             rhot(:) = trans(rho);
%             [imXt,imYt] = pol2cart(theta,rhot);
%             imt = interp2(imXt,imYt,im,imX,imY);
%             
%             ax.Children(i).CData = interp2(imXt,imYt,im,imX,imY);
%             
%             
%             imX = ax.Children(i).XData;
%             imY = ax.Children(i).YData';
%             [imX,imY] = meshgrid(ax.Children(i).XData,ax.Children(i).YData');
%             imXt = imX; imXt(:) = trans(imX);
%             imYt = imY; imYt(:) = trans(imY);
%             ax.Children(i).CData = interp2(imXt,imYt,im,imX,imY);
%             ax.Children(i).CData = interp2(trans(imX),trans(imY),im,imX,imY);
%             
%             Vq = interp2(X,Y,V,Xq,Yq)
%             
%             [theta,rho] = cart2pol(ax.Children(i).XData,ax.Children(i).YData');
%             [X,Y] = pol2cart(theta,trans(rho)');
%             ax.Children(i).XData = X;
%             ax.Children(i).YData = Y';
        otherwise
            error('X')
    end
end



a = 1;