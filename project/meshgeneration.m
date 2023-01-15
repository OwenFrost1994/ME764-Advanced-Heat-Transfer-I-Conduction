function [nElements,nNodes,connArray,coorNoUD,edgeArray,ElementsM,ElementsL]=meshgeneration(nElx_M,nElx_L,nEly,nElz,W,W1,H,T,order,nDof)
%this function works as the generator of mesh
%this function will calculate the number of nodes, number of elements, coordinates of nodes before the deformation,
%the connection array,in-matrix element index, in-layer element index
if(order==2)
    nNodes=(nElx_M+nElx_L+1)*(nEly+1)+(nElx_M+nElx_L+1)*(nEly)+(nElx_M+nElx_L)*(nEly+1);%number of nodes
    nElements=(nElx_M+nElx_L)*nEly;%number of elements
    elemW_M=(W-W1)/nElx_M;%width of element in matrix
    elemW_L=W1/nElx_L;%width of element in layer
    elemH=H/nEly;%hight of element
    %generate a matrix store the coordinates of nodes before the deformation
    coorNoUD=zeros(nNodes,2);
    for i=1:2*(nElx_M+nElx_L)+1
        if mod(i,2) ~= 0
            for j=1:2*nEly+1
                nid=(2*nEly+1)*(floor(i/2))+(nEly+1)*(floor(i/2))+j;
                if i<=nElx_M+1
                    coorNoUD(nid,1)=(i-1)*elemW_M/2;
                else if nElx_M+1<i && i<=nElx_M+2*nElx_L+1
                        coorNoUD(nid,1)=(W-W1)/2+(i-nElx_M-1)*elemW_L/2;
                    else
                        coorNoUD(nid,1)=(W+W1)/2+(i-nElx_M-2*nElx_L-1)*elemW_M/2;
                    end
                end
                coorNoUD(nid,2)=(2*nEly+1-j)*elemH/2;
            end
        else
            for j=1:nEly+1
                    nid=(2*nEly+1)*(floor(i/2))+(nEly+1)*(floor(i/2)-1)+j;
                    if i<=nElx_M+1
                        coorNoUD(nid,1)=(i-1)*elemW_M/2;
                    else if nElx_M+1<i && i<=nElx_M+2*nElx_L+1
                            coorNoUD(nid,1)=(W-W1)/2+(i-nElx_M-1)*elemW_L/2;
                        else
                            coorNoUD(nid,1)=(W+W1)/2+(i-nElx_M-2*nElx_L-1)*elemW_M/2;
                        end
                    end
                    coorNoUD(nid,2)=(nEly+1-j)*elemH;
            end
        end
    end
    % generate the connectivity array
    connArray=zeros(nElements,9);
    nEleM=1;
    nEleL=1;
    ElementsM=zeros(nElx_M*nEly,1);
    ElementsL=zeros(nElx_L*nEly,1);
    for i=1:nElx_M+nElx_L
        for j=1:nEly
            eid=j+(i-1)*nEly;
            connArray(eid,1)=1+2*j+(3*nEly+2)*(i-1);
            connArray(eid,2)=1+2*j+(3*nEly+2)*i;
            connArray(eid,3)=1+2*(j-1)+(3*nEly+2)*i;
            connArray(eid,4)=1+2*(j-1)+(3*nEly+2)*(i-1);
            connArray(eid,5)=1+j+(2*nEly+1)*i+(nEly+1)*(i-1);
            connArray(eid,6)=2*j+(3*nEly+2)*i;
            connArray(eid,7)=j+(2*nEly+1)*i+(nEly+1)*(i-1);
            connArray(eid,8)=2*j+(3*nEly+2)*(i-1);
            if nElx_M/2<i && i<=nElx_M/2+nElx_L
                ElementsL(nEleL,1)=eid;
                nEleL=nEleL+1;
                connArray(eid,9)=2;%in layer doman
            else
                ElementsM(nEleM,1)=eid;
                nEleM=nEleM+1;
                connArray(eid,9)=1;%in matirx doman
            end
        end
    end
    %the array store the node and boundary index of boundary edges
    %boundary 1 the left
    edgeArray=zeros(4*(nEly+nElx_M+nElx_L),3);
    for i=1:2*nEly
        edgeArray(i,1)=i;
        edgeArray(i,2)=i+1;
        edgeArray(i,3)=1;
    end
    %boundary 3 the right
    for i=1:2*nEly
        edgeArray(2*nEly+i,1)=(2*nEly+1+nEly+1)*(nElx_M+nElx_L)+i;
        edgeArray(2*nEly+i,2)=(2*nEly+1+nEly+1)*(nElx_M+nElx_L)+i+1;
        edgeArray(2*nEly+i,3)=3;
    end
    %boundary 2 the bottom
    for i=1:2*(nElx_M+nElx_L)
        if mod(i,2)==1
            edgeArray(4*nEly+i,1)=floor((i+1)/2)*(2*nEly+1)+floor((i-1)/2)*(nEly+1);
            edgeArray(4*nEly+i,2)=floor((i+1)/2)*(2*nEly+1)+floor((i-1)/2)*(nEly+1)+nEly+1;
            edgeArray(4*nEly+i,3)=2;
        else
            edgeArray(4*nEly+i,1)=floor((i)/2)*(3*nEly+2);
            edgeArray(4*nEly+i,2)=floor((i)/2)*(3*nEly+2)+2*nEly+1;
            edgeArray(4*nEly+i,3)=2;
        end
    end
    %boundary 4 the top
    for i=1:2*(nElx_M+nElx_L)
        if mod(i,2)==1
            edgeArray(4*nEly+2*(nElx_M+nElx_L)+i,1)=floor((i-1)/2)*(2*nEly+1)+floor((i-1)/2)*(nEly+1)+1;
            edgeArray(4*nEly+2*(nElx_M+nElx_L)+i,2)=floor((i-1)/2)*(2*nEly+1)+floor((i-1)/2)*(nEly+1)+2*nEly+2;
            edgeArray(4*nEly+2*(nElx_M+nElx_L)+i,3)=4;
        else
            edgeArray(4*nEly+2*(nElx_M+nElx_L)+i,1)=floor((i)/2)*(2*nEly+1)+(floor((i)/2)-1)*(nEly+1)+1;
            edgeArray(4*nEly+2*(nElx_M+nElx_L)+i,2)=floor((i)/2)*(2*nEly+1)+(floor((i)/2)-1)*(nEly+1)+1+nEly+1;
            edgeArray(4*nEly+2*(nElx_M+nElx_L)+i,3)=4;
        end
    end
end
if(order==1)
    nNodes=(nElx_M+nElx_L+1)*(nEly+1);%number of nodes
    nElements=(nElx_M+nElx_L)*nEly;%number of elements
    elemW_M=(W-W1)/nElx_M;%width of element in matrix
    elemH_M=H/nEly;%hight of element in matrix
    elemW_L=W1/nElx_L;%width of element in layer
    elemH=H/nEly;%hight of element in layer
    %generate a matrix store the coordinates of nodes before the deformation
    coorNoUD=zeros(nNodes,2);
    for i=1:nElx_M+nElx_L+1
        for j=1:nEly+1
            if i<=nElx_M/2+1
                coorNoUD(j+(i-1)*(nEly+1),1)=(i-1)*elemW_M;
            else if nElx_M/2+1<i && i<nElx_M/2+nElx_L+1
                    coorNoUD(j+(i-1)*(nEly+1),1)=(W-W1)/2+(i-nElx_M/2-1)*elemW_L;
                else
                    coorNoUD(j+(i-1)*(nEly+1),1)=(W+W1)/2+(i-nElx_M/2-nElx_L-1)*elemW_M;
                end
            end
            coorNoUD(j+(i-1)*(nEly+1),2)=(nEly+1-j)*elemH;
        end
    end
    % generate the connectivity array
    connArray=zeros(nElements,5);
    nEleM=1;
    nEleL=1;
    ElementsM=zeros(nElx_M*nEly,1);
    ElementsL=zeros(nElx_L*nEly,1);
    for i=1:nElx_M+nElx_L
        for j=1:nEly
            eid=j+(i-1)*nEly;
            connArray(eid,1)=1+j+(i-1)*(nEly+1);
            connArray(eid,2)=1+j+i*(nEly+1);
            connArray(eid,3)=j+i*(nEly+1);
            connArray(eid,4)=j+(i-1)*(nEly+1);
            if nElx_M/2<i && i<=nElx_M/2+nElx_L
                ElementsL(nEleL,1)=eid;
                nEleL=nEleL+1;
                connArray(eid,5)=2;%in layer doman
            else
                ElementsM(nEleM,1)=eid;
                nEleM=nEleM+1;
                connArray(eid,5)=1;%in matrix doman
            end
        end
    end
    
    %the array store the node and boundary index of boundary edges
    %boundary 1 the left
    edgeArray=zeros(2*(nEly+nElx_M+nElx_L),3);
    for i=1:nEly
        edgeArray(i,1)=i;
        edgeArray(i,2)=i+1;
        edgeArray(i,3)=1;
    end
    %boundary 3 the right
    for i=1:nEly
        edgeArray(nEly+i,1)=(nEly+1)*(nElx_M+nElx_L)+i;
        edgeArray(nEly+i,2)=(nEly+1)*(nElx_M+nElx_L)+i+1;
        edgeArray(nEly+i,3)=3;
    end
    %boundary 2 the bottom
    for i=1:nElx_M+nElx_L
        edgeArray(2*nEly+i,1)=i*(nEly+1);
        edgeArray(2*nEly+i,2)=(i+1)*(nEly+1);
        edgeArray(2*nEly+i,3)=2;
    end
    %boundary 4 the top
    for i=1:nElx_M+nElx_L
        edgeArray(2*nEly+nElx_M+nElx_L+i,1)=(i-1)*(nEly+1)+1;
        edgeArray(2*nEly+nElx_M+nElx_L+i,2)=i*(nEly+1)+1;
        edgeArray(2*nEly+nElx_M+nElx_L+i,3)=4;
    end
end
end