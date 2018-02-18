function boundaries = plotBoundaries(data)

% Boundary
fieldNames = fieldnames(data.PGDmeshes(1).elementFaceInfo);
nOfBoundaries = length(fieldNames);
boundaries = cell(1,nOfBoundaries);
for i=1:nOfBoundaries
    iboundary = fieldNames{i};
    nOfElements = size(data.PGDmeshes(1).elementFaceInfo.(iboundary),1);
    boundaries{i} = zeros(nOfElements,1);
    hold on
    for j=1:nOfElements
        Tb = data.PGDmeshes(1).elementFaceInfo.(iboundary)(j,:);
        Tf = data.PGDmeshes(1).T.all(Tb(1),data.PGDmeshes(1).referenceElement.faceNodes(Tb(2),:));
        Xf = data.PGDmeshes(1).X(Tf,:);
        boundaries{i}(j) = plot(Xf(:,1),Xf(:,2),'-r','LineWidth',3,'Visible','off');
    end
%     boundaries{i} = plot(Xf(:,1),Xf(:,2),'-r','LineWidth',10,'Marker','.','MarkerSize',10,'Visible','off');
    hold off
end
axis tight