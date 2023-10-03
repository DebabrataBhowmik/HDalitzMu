#include <vector>
#include <algorithm>

std::vector<float> sigmaEff(std::vector<float> v, const float threshold){
    std::sort(v.begin(), v.end());

    int total = v.size();
    int max = (int)(threshold * total);

    std::vector<float>  start;
    std::vector<float>  stop;
    std::vector<float>  width;

    int i = 0;
    while (i != total-1){
        int count = 0;
        int j = i;
        while (j != total-1 && count < max){
            ++count;
            ++j;
        }
        if (j != total-1){
            start.push_back(v[i]);
            stop .push_back(v[j]);
            width.push_back(v[j] - v[i]);
        }
        ++i;
    }

    float minwidth = *min_element(width.begin(), width.end()) * 0.5;
    auto pos = min_element(width.begin(), width.end()) - width.begin();
    float xmin = start[pos];
    float xmax = stop[pos];

    std::vector<float> sigma = {xmin, xmax, minwidth};
    return sigma;
}