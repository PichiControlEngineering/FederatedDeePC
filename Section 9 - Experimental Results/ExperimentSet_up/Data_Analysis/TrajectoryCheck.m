dataPath = convertCharsToStrings(findFolder("Experiment Data"));
Pato_List = [1, 7, 11, 12, 14, 17];
y_list = zeros(24001, length(Pato_List));

for i = 1:length(Pato_List)
    Pato_i = Pato_List(i);
    load(dataPath +  "\"+ Pato_i + "\PatoResp2_" + 250 + "Hz.mat", 'data')
    y_i = data{2}.Values.Data;
    y_list(:,i) = y_i;
end

figure()
plot(y_list(1:300,[1,2,3,4,5])); legend({"1", "2", "3", "4", "5"}); hold on
plot(mean(y_list(1:300,:),2), 'k--')
