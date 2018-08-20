function [G,n,k,H]=rref_mod2()

H_14_7;

fileName = ['' num2str(size(H,2)) '_' num2str(size(H,1)) '.bin'];

fileID = fopen(fileName, 'w');

fwrite(fileID, size(H,2), 'integer*4');
fwrite(fileID, size(H,1), 'integer*4');

for i = 1:size(H,1)
    
    for j = 1:size(H,2)
        fwrite(fileID, H(i,j));
    end    
end

n = min(size(H));

% Swap rows so that a diagonal of ones is formed
for i = 1:(n-1)
    % Check if the (i,i) element is a 1, if it's not switch with the next row
    % that has a leading 1 in the (i,i)th position
    if H(i,i) == 0
        % Check for the next row with a zero
        for j = (i+1):size(H,1)
            if H(j,i) == 1
                H([i j],:) = H([j i],:);
                break;
            end % if
        end % for
    end % if
    
    % clear any ones below the leading on in the rest of the column
    for j = (i+1):size(H,1)
        if H(j,i) == 1
            H(j,:) = mod(H(j,:) + H(i,:),2);
        end
    end
    
end % for

 % Remove any rows of all zero
 H(all(H==0,2),:) = [];

% % Swap any columns so that the identity matrix is on the left
% for i = 1:n
%     % If we don't have the identity matrix yet
%     if(i > size(H,1))
%         break;
%     end
%     
%     if H(i,i) == 0
%         %Find the column to swap it with 
%         for j = (i+1):size(H,2)
%             if H(i,j) == 1
%                 H(:,[i j]) = H(:,[j i]);
%                 swap = 1;
%                 break;
%             end
%         end
%     end
%     
%     % clear any ones below the leading one in the rest of the column
%     if i ~= n
%         for j = (i+1):size(H,1)
%             if H(j,i) == 1
%                 H(j,:) = mod(H(j,:) + H(i,:),2);
%             end
%         end
%     end
%     
%     % Remove any rows of all zero
%     H(all(H==0,2),:) = [];
% end

% clear the columns after the leading one's
for i = 2:size(H,1)
    
    for j = 1:i-1
        if H(j,i) == 1
            H(j,:) = mod(H(j,:) + H(i,:),2);
        end
    end
end

% Remove any rows of all zero
H(all(H==0,2),:) = []

G = gen2par(H)

[k, n] = size(G);

fwrite(fileID,size(G,1),'integer*4');

for i = 1:size(G,1)
    
    for j = 1:size(G,2)
        fwrite(fileID, G(i,j));
    end
    
end

fwrite(fileID,k,'integer*4');

fclose(fileID);



