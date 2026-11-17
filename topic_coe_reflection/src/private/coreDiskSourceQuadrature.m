function [positions,weights]=coreDiskSourceQuadrature(diameter,radialNodes,azimuthNodes)
%coreDiskSourceQuadrature 均匀非相干圆盘的Gauss面积积分，t=r²/R²∈[0,1]。
    validateattributes(diameter,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(radialNodes,{'numeric'},{'scalar','integer','positive'});
    validateattributes(azimuthNodes,{'numeric'},{'scalar','integer','>=',4});
    if diameter==0, positions=[0,0]; weights=1; return; end
    if radialNodes==1
        nodes=0.5; radialWeights=1;
    else
        j=(1:radialNodes-1).'; beta=j./sqrt(4*j.^2-1);
        [v,d]=eig(diag(beta,1)+diag(beta,-1));
        [nodes,order]=sort((diag(d)+1)/2); radialWeights=v(1,order).^2;
    end
    phi=2*pi*(0:azimuthNodes-1)/azimuthNodes;
    radius=diameter/2*sqrt(nodes(:));
    sx=radius*cos(phi); sy=radius*sin(phi);
    positions=[sx(:),sy(:)];
    weightMap=repmat(radialWeights(:)/azimuthNodes,1,azimuthNodes);
    weights=reshape(weightMap,1,[]);
    weights=weights/sum(weights);
end
