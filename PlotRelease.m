figure;
 
labels = {"Regular", "Mature", "Immature 1", "Immature 2", "Immature 3", "Immature 4", "Immature 5"};

cmap = jet(7);

times = 0:1e-4:1;

for i = 1:length(labels)
    plot(data{i}(1,:) * 1000, data{i}(2, :), "Color", cmap(i,:))
    hold on
end

legend(labels, "Location", "northwest");
xlabel("Peak Voltage (mV)");
ylabel("Vesicle Release (Hz)");
title("Vesicle Release Per Second");

figure;
for i = 1:length(labels)
    plot(times(2:end), cumsum(cumu_release{i}))
    hold on 
end

legend(labels, "Location", "northwest");
xlabel("Time (s)");
ylabel("Vesicle Release")
title("Cumulative Vesicle Release at -25 mV");