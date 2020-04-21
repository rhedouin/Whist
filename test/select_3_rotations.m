%
clear
close all

load('BrainSample2_rotations_ref_2_orientations.mat')
load('20_fiber_orientations_rotations.mat')

C = nchoosek(1:9,3)

for k = 1:size(C,1)
% for k = [66, 47]

    for m = 1:20
        for l=1:3
            display(['k: ' num2str(k) ', l: ' num2str(l) ', m: ' num2str(m)])
            current_rotation = rotations(:,:,C(k,l));
            fiber_rotated = current_rotation*fiber_directions(m,:)';
            theta(l,m) = acos(abs(fiber_rotated(3)));
            all_theta(k,l,m) = theta(l,m);
        end
        theta_sorted(:,m) = sort(theta(:,m));
        diff_theta(m) = theta_sorted(end,m) -  theta_sorted(1,m);
        
    end
%     
%         figure
%         plot(theta_sorted)
        
    avg_diff_theta(k) = mean(diff_theta);
end         
[a, b] = max(avg_diff_theta)
best_rotation = C(b,:)
keyboard;

        


