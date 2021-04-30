function TBL = manovaX(x)

withinFac1 = repmat({'real' 'imag'}',[size(x,2) 1]);
withinFac2 = repmat({'grat1' 'grat2' 'plaid'},[2 1]); withinFac2 = withinFac2(:);
withinDesign = table(withinFac1, withinFac2);
withinDesign.Properties.VariableNames = {'complex' 'stimCond'};
withinModel = 'stimCond*complex';

t = table(real(x(:,1)),imag(x(:,1)),real(x(:,2)),imag(x(:,2)),real(x(:,3)),imag(x(:,3)));
t.Properties.VariableNames = {'grat1Real' 'grat1Imag' 'grat2Real' 'grat2Imag' 'plaidReal' 'plaidImag'};

rm = fitrm(t,'grat1Real-plaidImag~1','WithinDesign',withinDesign,'WithinModel',withinModel);
[TBL,A,C_MV,D,withinNames,betweenNames] = manova2(rm,withinModel);