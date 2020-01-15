function cverror = spdaux_smoothcvsingle(input, x, hval, inds)

% given index 'inds', compute cross-validation error

uind = unique(inds);
K = length(uind);
p = input.size(1);
cverror = 0;
for i=1:K
    idtest  = find(inds==uind(i));
    idtrain = find(inds~=uind(i));
    
    Ytrain = struct();
    Ytrain.name = 'spd';
    Ytrain.size = [p,p,length(idtrain)];
    Ytrain.data = input.data(:,:,idtrain);
    xtrain = x(idtrain);
    xtest  = x(idtest);
    
    % let's try to predict
    yeval  = spdaux_smootheval(Ytrain, hval, xtrain);
    Ypred  = zeros(p,p,length(idtest));
    for j=1:length(idtest)
        Ypred(:,:,j) = yeval(xtest(j));
    end
    Ytest  = input.data(:,:,idtest);
    
    % compute added cverror
    cvadded = 0;
    for j=1:length(idtest)
        cvadded = cvadded + norm(Ypred(:,:,j)-Ytest(:,:,j),'fro');
    end
    cvadded = cvadded/length(idtest);
    
    % update
    cverror = cverror + cvadded;
end
end
