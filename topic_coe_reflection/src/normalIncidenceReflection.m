function coefficient = normalIncidenceReflection(incidentIndex, refractiveIndex, extinctionCoefficient)
%normalIncidenceReflection 法向入射、单界面、半无限材料的复振幅 Fresnel 反射系数。
% 采用 N=n+i*k（时间因子 exp(-i*omega*t)）；不包含薄膜干涉或高度产生的传播相位。
    arguments
        incidentIndex (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
        refractiveIndex double {mustBeReal, mustBeFinite, mustBePositive}
        extinctionCoefficient double {mustBeReal, mustBeFinite, mustBeNonnegative}
    end
    assert(isequal(size(refractiveIndex), size(extinctionCoefficient)), ...
        'COEReflection:IndexSize', 'n 和 k 必须具有相同尺寸。');
    materialIndex = refractiveIndex + 1i*extinctionCoefficient;
    coefficient = (incidentIndex - materialIndex)./(incidentIndex + materialIndex);
end
