clear
close all

base_folder = '/project/3015069.04/signal_components/multi_orientations/BrainSample2/';

old = load([base_folder 'without_dispersion_save/FVF60_N400_train1/Signal_FVF60_train1.mat'])
new = load([base_folder 'without_dispersion/FVF60_N400_train1/Signal_FVF60_train1.mat'])

for k = [2 15]
    for l = [4 8]
        old_signal = old.SignalComponent(1, 3, 10, l);
        new_signal = new.SignalComponent(1, 1, 1, 4, l);
        
        old_dir = old_signal.current_dir
        new_dir = new_signal.current_dir
        
        time = old.time;
        
        figure
        subplot(221)
        plot(time, old_signal.Axon)
        hold on
        plot(time, new_signal.Axon)
        ylim([-0.2 0.5])
        
        subplot(222)
        plot(time, old_signal.Myelin)
        hold on
        plot(time, new_signal.Myelin)
              ylim([-0.2 0.5])
  
        subplot(223)
        plot(time, old_signal.Extra)
        hold on
        plot(time, new_signal.Extra)
        ylim([-0.2 0.5])

    end
end

% for k = [2 15]
%     for l = [4 8]
%         old_signal = old.SignalComponent(1, 3, k + 3, l);
%         new_signal = new.SignalComponent(1, 1, 1, k, l);
%         
%         old_dir = old_signal.current_dir
%         new_dir = new_signal.current_dir
%         
%         time = old.time;
%         
%         figure
%         subplot(221)
%         plot(time, old_signal.Axon)
%         hold on
%         plot(time, new_signal.Axon)
%         
%         subplot(222)
%         plot(time, old_signal.Myelin)
%         hold on
%         plot(time, new_signal.Myelin)
%         
%         subplot(223)
%         plot(time, old_signal.Extra)
%         hold on
%         plot(time, new_signal.Extra)
%     end
% end