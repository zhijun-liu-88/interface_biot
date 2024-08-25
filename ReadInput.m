
function [] = ReadInput(InputFile)

% global JointBoundLine NumHolePoint CoorHolePoint MeshShape
% global SizeETri HSizeRTri VSizeRTri MinXCenterRTriMC MinYCenterRTriMC
% global HSizeQuad VSizeQuad MinXCenterQuadMC MinYCenterQuadMC
% 
% InputID = fopen(InputFile, 'r');
% 
% %fprintf("InputID =%d\n", InputID);
% 
% TempStr = fgetl(InputID);
% MeshShape = fgetl(InputID);
% 
% TempStr = fgetl(InputID);
% NumBoundSegs = fscanf(InputID, '%d');
% 
% TempStr = fgetl(InputID);
% formatSpec = '%f %f %f %f';
% sizeA = [4  NumBoundSegs];
% JointBoundLine = fscanf(InputID, formatSpec, sizeA);
% JointBoundLine = JointBoundLine.';
% TempStr = fgetl(InputID);
% 
% TempStr = fgetl(InputID);
% NumHolePoint = fscanf(InputID, '%d');
% 
% if NumHolePoint > 0
%     TempStr = fgetl(InputID);
%     formatSpec = '%f %f';
%     sizeA = [2 NumHolePoint];
%     CoorHolePoint = fscanf(InputID, formatSpec, sizeA);
%     CoorHolePoint = CoorHolePoint.';
%     TempStr = fgetl(InputID);
% end 
% 
% TempStr = fgetl(InputID);
% formatSpec = '%f';
% SizeETri = fscanf(InputID, formatSpec);
% 
% TempStr = fgetl(InputID);
% formatSpec = '%f %f';
% sizeA = [2 1];
% TempLine = fscanf(InputID, formatSpec, sizeA);
% HSizeRTri = TempLine(1); 
% VSizeRTri = TempLine(2);
% TempStr = fgetl(InputID);
% 
% TempStr = fgetl(InputID);
% formatSpec = '%f %f';
% sizeA = [2 1];
% TempLine = fscanf(InputID, formatSpec, sizeA);
% MinXCenterRTriMC = TempLine(1); 
% MinYCenterRTriMC = TempLine(2);
% TempStr = fgetl(InputID);
% 
% TempStr = fgetl(InputID);
% formatSpec = '%f %f';
% sizeA = [2 1];
% TempLine = fscanf(InputID, formatSpec, sizeA);
% MinXCenterQuadMC = TempLine(1); 
% MinYCenterQuadMC = TempLine(2);
% TempStr = fgetl(InputID);
% 
% TempStr = fgetl(InputID);
% formatSpec = '%f %f';
% sizeA = [2 1];
% TempLine = fscanf(InputID, formatSpec, sizeA);
% HSizeQuad = TempLine(1); 
% VSizeQuad = TempLine(2);
% TempStr = fgetl(InputID);
% 
% fclose(InputID);





global JointBoundLine NumHolePoint CoorHolePoint MeshShape
global SizeETri HSizeRTri VSizeRTri MinXCenterRTriMC MinYCenterRTriMC
global HSizeQuad VSizeQuad MinXCenterQuadMC MinYCenterQuadMC

InputID = fopen(InputFile, 'r');

%fprintf("InputID =%d\n", InputID);

TempStr = fgetl(InputID);
MeshShape = fgetl(InputID);

TempStr = fgetl(InputID);
NumBoundSegs = fscanf(InputID, '%d');

TempStr = fgetl(InputID);
formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f';
sizeA = [4  NumBoundSegs];
JointBoundLine = fscanf(InputID, formatSpec, sizeA);
JointBoundLine = JointBoundLine.';
TempStr = fgetl(InputID);

TempStr = fgetl(InputID);
NumHolePoint = fscanf(InputID, '%d');

if NumHolePoint > 0
    TempStr = fgetl(InputID);
    formatSpec = '%f %f';
    sizeA = [2 NumHolePoint];
    CoorHolePoint = fscanf(InputID, formatSpec, sizeA);
    CoorHolePoint = CoorHolePoint.';
    TempStr = fgetl(InputID);
end 

TempStr = fgetl(InputID);
formatSpec = '%f';
SizeETri = fscanf(InputID, formatSpec);

TempStr = fgetl(InputID);
formatSpec = '%f %f';
sizeA = [2 1];
TempLine = fscanf(InputID, formatSpec, sizeA);
HSizeRTri = TempLine(1); 
VSizeRTri = TempLine(2);
TempStr = fgetl(InputID);

TempStr = fgetl(InputID);
formatSpec = '%f %f';
sizeA = [2 1];
TempLine = fscanf(InputID, formatSpec, sizeA);
MinXCenterRTriMC = TempLine(1); 
MinYCenterRTriMC = TempLine(2);
TempStr = fgetl(InputID);

TempStr = fgetl(InputID);
formatSpec = '%f %f';
sizeA = [2 1];
TempLine = fscanf(InputID, formatSpec, sizeA);
MinXCenterQuadMC = TempLine(1); 
MinYCenterQuadMC = TempLine(2);
TempStr = fgetl(InputID);

TempStr = fgetl(InputID);
formatSpec = '%f %f';
sizeA = [2 1];
TempLine = fscanf(InputID, formatSpec, sizeA);
HSizeQuad = TempLine(1); 
VSizeQuad = TempLine(2);
TempStr = fgetl(InputID);



% HSizeQuad = 10 / 3;
% VSizeQuad = 10 / 3;



fclose(InputID);





