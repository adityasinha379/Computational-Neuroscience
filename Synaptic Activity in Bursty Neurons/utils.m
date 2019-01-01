classdef utils
    methods(Static)
        function [spike_states] = poissonSpikes(t, rate, nTrials)
            if nargin < 3
                nTrials = 1;
            end
            dt = diff(t(1:2)); % unit in s
            spike_states = rand([nTrials, length(t)]) < rate*dt;
        end
       
        function plotRaster(t, spike_state,cs)
            spike_state = logical(spike_state);
            nTrials = size(spike_state,1);
            if nargin<3
                cs = zeros(nTrials,3);
            end
            hold all;
            for tr = 1:nTrials
                tk = t(spike_state(tr,:)); tk = tk(:);
                plot([tk,tk]',[ones(size(tk))*tr-0.4,ones(size(tk))*tr+0.4]','k');
                plot([min(t),min(t)], [tr-0.5,tr+0.5],'Color',cs(tr,:),'LineWidth',3);
            end
            ylim([0,size(spike_state, 1)+1]);
            xlim([min(t),max(t)]);
        end
        
        function plotDynamic(t,Iext,V,R,isVideo)
            if nargin<5
                isVideo=0;
            end
            ax1 = subplot(221); hold on;
            title('$I_{ext}$');
            plot(t,Iext);
            ax2 = subplot(223); hold on;
            title('V');
            plot(t,V);
            ax3 = subplot(2,2,[2,4]); hold on;
            title('V vs. R');
            plot(V,R);
            
            if isVideo
                for tt = 1:100:length(t)
                    l1 = plot(ax1,t(tt),Iext(tt),'ko');
                    set(l1,'markerfacecolor','g');
                    l2 = plot(ax2,t(tt),V(tt),'ko');
                    set(l2,'markerfacecolor','g');
                    l3 = plot(ax3,V(tt),R(tt),'ko');
                    set(l3,'markerfacecolor','g');
                    pause(0.01);
                    delete([l1,l2,l3]);
                end
            end
        end
    end
end