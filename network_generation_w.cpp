#include <iostream>
#include <fstream>
#include <locale>
#include <codecvt>
#include <map>
#include <vector>
#include <algorithm>
#include <cassert>
#include <set>
#include <queue>
#include <deque>
#include <thread>
#include <atomic>
#include <random>
using namespace std;

static const std::wstring punc = L"，《》。、「！？」：；-’‘“”…■『』〔〕（）．:4";
/*
inline bool not_punc(const std::wstring &str1) {
    if (str1.find_first_of(punc) != static_cast<size_t>(-1)) {
        return false;
    }
    return true;
}

inline bool not_punc(wchar_t c1) {
    if (punc.find_first_of(c1) != static_cast<size_t>(-1)) {
        return false;
    }
    return true;
}
*/
template<typename T>
inline bool not_punc(const T &a) {
    if (punc.find_first_of(a) != static_cast<size_t>(-1)) {
        return false;
    }
    return true;
}

int main (int argc, char *argv[]) {
    if (argc <= 1) {
        cerr << "Provide File Name" << endl;
        return 0;
    }

    const locale utf8_locale(locale(), new codecvt_utf8<wchar_t>());
    wcout.imbue(utf8_locale);
    std::wstring str;

    wifstream word_ifs("words.txt", wifstream::in);
    word_ifs.imbue(utf8_locale);

    size_t count = 0;
    map<std::wstring,size_t> dictionary;
    map<std::wstring,size_t> id;
    vector<std::wstring> r_id;
    while (word_ifs >> str) {
        r_id.push_back(str);
        id[str] = count++;
        dictionary[str] = 0;
        word_ifs >> str;
        word_ifs >> str;
    }

    vector<map<std::wstring,size_t>> E;
    E.resize(count);

    set<pair<size_t,size_t>> edge;
    wifstream ifs (argv[1], wifstream::in);
    ifs.imbue(utf8_locale);

    while (ifs >> str) {
        size_t start_1 = 0;
        size_t len_1 = 1;
        while (dictionary.find(str.substr(start_1,++len_1)) != dictionary.end() && start_1 + len_1 <= str.length()) {}

        std::wstring str1 = str.substr(start_1,--len_1);
        while (start_1 + len_1 < str.length()) {
            //wcout << str1 << endl;
            size_t start_2 = start_1 + len_1;
            size_t len_2 = 1;
            while (dictionary.find(str.substr(start_2,++len_2)) != dictionary.end() && start_2 + len_2 <= str.length()) {}

            std::wstring str2 = str.substr(start_2,--len_2);
            if (not_punc(str1) && not_punc(str2) && str1 != str2 && dictionary.find(str1) != dictionary.end() && dictionary.find(str2) != dictionary.end()) {
                ++E[id[str1]][str2];
                ++E[id[str2]][str1];
                if (edge.find(make_pair(id[str2],id[str1])) == edge.end())
                    edge.insert(make_pair(id[str1],id[str2]));
            }

            start_1 = start_2;
            len_1 = len_2;
            str1 = std::move(str2);
        }
    }

    // Find most powerful node
    vector<pair<size_t,size_t>> Esize;
    for (size_t i = 0 ; i < E.size() ; ++i) {
        Esize.push_back(make_pair(E[i].size(),i));
    }
    sort(Esize.begin(),Esize.end());

    // Find the largest component
    vector<bool> lc;
    lc.resize(E.size());
    queue<size_t> q;
    q.push(Esize[E.size() - 1].second);
    while (!q.empty()) {
        size_t v = q.front();
        q.pop();
        if (lc[v]) {
            continue;
        }
        lc[v] = true;
        for (auto e : E[v]) {
            q.push(id[e.first]);
        }
    }
    map<std::wstring,size_t> id_c;
    vector<std::wstring> r_id_c;
    vector<map<std::wstring,size_t>> E_c;
    for (size_t i = 0 ; i < E.size() ; ++i) {
        if (lc[i]) {
            id_c[r_id[i]] = r_id_c.size();
            r_id_c.push_back(r_id[i]);
            map<std::wstring,size_t> temp;
            for (auto e : E[i]) {
                temp[e.first] = e.second;
            }
            E_c.push_back(temp);
        }
    }
    id = id_c;
    r_id = r_id_c;
    E = E_c;
    edge.clear();
    for (size_t i = 0 ; i < E.size() ; ++i) {
        for (auto e : E[i]) {
            if (edge.find(make_pair(id[e.first],i)) == edge.end()) {
                edge.insert(make_pair(i,id[e.first]));
            }
        }
    }
    Esize.clear();
    for (size_t i = 0 ; i < E.size() ; ++i) {
        Esize.push_back(make_pair(E[i].size(),i));
    }
    sort(Esize.begin(),Esize.end());

/*
    for (size_t i = 0 ; i < r_id.size() ; ++i) {
        // Node List
        wcout << i << L"," << r_id[i] << endl;
    }
*/
/*
    for (auto item : id) {
        // Node List
        wcout << item.first << L"\t" << item.second << endl;
    }
*/
/*
    count = 0;
    for (size_t i = 0 ; i < E.size() ; ++i) {
        // Edge List
        for (auto p : E[i]) {
            if (i < id[p.first]) {
                wcout << count++ << L"," << i << L"," << id[p.first] << L"," << p.second << endl;
            }
        }
    }
*/
/*
    for (auto e : Esize) {
        wcout << r_id[e.second] << endl;
    }
*/

    // Data Processing
    // Calculate total number of links
    size_t total_number_link = 0;
    for (auto v : E) {
        total_number_link += v.size();
    }
    assert((total_number_link & 1) == 0);
    total_number_link >>= 1;
/*
    {
        // degree correlation measured by knn(k)
        map<size_t,double> knn;
        size_t start = 0;
        double cml = 0;
        for (size_t i = 0 ; i < Esize.size() ; ++i) {
            size_t temp = 0;
            for (auto v : E[Esize[i].second]) {
                temp += E[id[v.first]].size();
            }
            cml += static_cast<double>(temp) / static_cast<double>(E[Esize[i].second].size());

            if (i + 1 >= Esize.size() || Esize[i + 1].first != Esize[start].first) {
                knn[Esize[start].first] = cml / static_cast<double>(i - start + 1);
                start = i + 1;
                cml = 0;
            }
        }
        for (auto v : knn) {
            wcout << L"(" << v.first << L"," << v.second << L")" << endl;
        }
    }
*/
/*
    {
        // Assortative coefficient alpha
        size_t sum_x = 0;
        size_t sum_xy = 0;
        size_t sum_x_2 = 0;
        for (auto v : E) {
            size_t temp = 0;
            for (auto u : v) {
                temp += E[id[u.first]].size();
                sum_x += v.size();
                sum_x_2 += v.size() * v.size();
            }
            sum_xy += temp * v.size();
        }
        double alpha = (static_cast<double>(total_number_link * static_cast<size_t>(2) * sum_xy) - static_cast<double>(sum_x * sum_x)) / (static_cast<double>(total_number_link * static_cast<size_t>(2) * sum_x_2) - static_cast<double>(sum_x * sum_x));
        wcout << L"alpha: " << alpha << endl;
    }


    {
        // Another method for alpha
        double sum_kiki = 0;
        double sum_ki = 0;
        double sum_ki_2 = 0;
        for (auto e : edge) {
            sum_kiki += E[e.first].size() * E[e.second].size();
            sum_ki += E[e.first].size() + E[e.second].size();
            sum_ki_2 +=  E[e.first].size() *  E[e.first].size() + E[e.second].size() * E[e.second].size();
        }
        double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        wcout << L"alpha: " << alpha << endl;
    }
*/
/*
    {
        // Calculate Rich club coefficient
        size_t cum_E = 0;
        vector<size_t> nodes_account = {Esize[Esize.size() - 1].second};
        for (int i = Esize.size() - 2 ; i >= 0 ; --i) {
            size_t u = Esize[i].second;
            if (E[u].size() == 0) {
                break;
            }
            for (auto v : nodes_account) {
                if (edge.find(make_pair(u,v)) != edge.end() || edge.find(make_pair(v,u)) != edge.end()) {
                    ++cum_E;
                }
            }
            nodes_account.push_back(u);
            double temp = cum_E;
            temp *= 2.0L / static_cast<double>(nodes_account.size() * (nodes_account.size() - 1));
            wcout << L"(" << nodes_account.size() << L"," << temp << L")" << endl;
        }
    }
*/
/*
    {
        // Link Rewiring (MAX assortative)
        vector<pair<size_t,size_t>> edge_copy;
        for (auto e : edge) {
            edge_copy.push_back(e);
        }
        for (auto it1 = edge_copy.begin() ; it1 != edge_copy.end() ; ++it1) {
            auto it2 = it1;
            for (++it2 ; it2 != edge_copy.end() ; ++it2) {
                vector<pair<size_t,size_t>> temp;
                temp.push_back(make_pair(E[it1->first].size(),it1->first));
                temp.push_back(make_pair(E[it1->second].size(),it1->second));
                temp.push_back(make_pair(E[it2->first].size(),it2->first));
                temp.push_back(make_pair(E[it2->second].size(),it2->second));
                sort(temp.begin(),temp.end());
                if (temp[0].second != temp[1].second && temp[2].second != temp[3].second) {
                    it1->first = temp[0].second;
                    it1->second = temp[1].second;
                    it2->first = temp[2].second;
                    it2->second = temp[3].second;
                }
            }
        }
        double sum_kiki = 0;
        double sum_ki = 0;
        double sum_ki_2 = 0;
        for (auto e : edge_copy) {
            sum_kiki += E[e.first].size() * E[e.second].size();
            sum_ki += E[e.first].size() + E[e.second].size();
            sum_ki_2 +=  E[e.first].size() *  E[e.first].size() + E[e.second].size() * E[e.second].size();
        }
        double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        wcout << L"alpha: " << alpha << endl;
    }

    {
        // Link Rewiring (MAX disassortative)
        vector<pair<size_t,size_t>> edge_copy;
        for (auto e : edge) {
            edge_copy.push_back(e);
        }
        for (auto it1 = edge_copy.begin() ; it1 != edge_copy.end() ; ++it1) {
            auto it2 = it1;
            for (++it2 ; it2 != edge_copy.end() ; ++it2) {
                vector<pair<size_t,size_t>> temp;
                temp.push_back(make_pair(E[it1->first].size(),it1->first));
                temp.push_back(make_pair(E[it1->second].size(),it1->second));
                temp.push_back(make_pair(E[it2->first].size(),it2->first));
                temp.push_back(make_pair(E[it2->second].size(),it2->second));
                sort(temp.begin(),temp.end());
                if (temp[0].second != temp[3].second && temp[1].second != temp[2].second) {
                    it1->first = temp[0].second;
                    it1->second = temp[3].second;
                    it2->first = temp[1].second;
                    it2->second = temp[2].second;
                }
            }
        }
        double sum_kiki = 0;
        double sum_ki = 0;
        double sum_ki_2 = 0;
        for (auto e : edge_copy) {
            sum_kiki += E[e.first].size() * E[e.second].size();
            sum_ki += E[e.first].size() + E[e.second].size();
            sum_ki_2 +=  E[e.first].size() *  E[e.first].size() + E[e.second].size() * E[e.second].size();
        }
        double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        wcout << L"alpha: " << alpha << endl;
    }

    {
        // Link Rewiring (MAX random)
        const size_t N = 30;
        default_random_engine generator(time(NULL));
        uniform_int_distribution<int> distribution(1,3);
        double alpha_cum = 0;
        for (size_t rr = 0 ; rr < N ; ++rr) {
            vector<pair<size_t,size_t>> edge_copy;
            for (auto e : edge) {
                edge_copy.push_back(e);
            }
            for (auto it1 = edge_copy.begin() ; it1 != edge_copy.end() ; ++it1) {
                auto it2 = it1;
                for (++it2 ; it2 != edge_copy.end() ; ++it2) {
                    vector<pair<size_t,size_t>> temp;
                    temp.push_back(make_pair(E[it1->first].size(),it1->first));
                    temp.push_back(make_pair(E[it1->second].size(),it1->second));
                    temp.push_back(make_pair(E[it2->first].size(),it2->first));
                    temp.push_back(make_pair(E[it2->second].size(),it2->second));
                    std::swap(temp[1],temp[distribution(generator)]);
                    if (temp[0].second != temp[1].second && temp[2].second != temp[3].second) {
                        it1->first = temp[0].second;
                        it1->second = temp[1].second;
                        it2->first = temp[2].second;
                        it2->second = temp[3].second;
                    }
                }
            }
            double sum_kiki = 0;
            double sum_ki = 0;
            double sum_ki_2 = 0;
            for (auto e : edge_copy) {
                sum_kiki += E[e.first].size() * E[e.second].size();
                sum_ki += E[e.first].size() + E[e.second].size();
                sum_ki_2 +=  E[e.first].size() *  E[e.first].size() + E[e.second].size() * E[e.second].size();
            }
            double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
            alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
            wcout << L"R alpha: " << alpha << endl;
            alpha_cum += alpha;
        }
        wcout << L"Avg alpha: " << alpha_cum / static_cast<double>(N) << endl;
    }
*/

/*    
    {
        // Second order average assortative coef
        vector<double> K;
        K.resize(E.size());
        for (size_t i = 0 ; i < E.size() ; ++i) {
            K[i] = 0;
            for (auto u : E[i]) {
                K[i] += E[id[u.first]].size();
            }
            K[i] /= static_cast<double>(E[i].size());
        }
        double sum_kiki = 0;
        double sum_ki = 0;
        double sum_ki_2 = 0;
        for (auto e : edge) {
            sum_kiki += K[e.first] * K[e.second];
            sum_ki += K[e.first] + K[e.second];
            sum_ki_2 += K[e.first] * K[e.first] + K[e.second] * K[e.second];
        }
        double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        wcout << L"Second order alpha: " << alpha << endl;
    }


    {
        // Link Rewiring (second order MAX assortative)
        vector<double> K;
        K.resize(E.size());
        for (size_t i = 0 ; i < E.size() ; ++i) {
            K[i] = 0;
            for (auto u : E[i]) {
                K[i] += E[id[u.first]].size();
            }
            K[i] /= static_cast<double>(E[i].size());
        }
        vector<pair<size_t,size_t>> edge_copy;
        for (auto e : edge) {
            edge_copy.push_back(e);
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto it1 = edge_copy.begin() ; it1 != edge_copy.end() ; ++it1) {
                auto it2 = it1;
                for (++it2 ; it2 != edge_copy.end() ; ++it2) {
                    if (E[it1->first].size() == E[it2->first].size()) {
                        if ((K[it1->second] < K[it2->second] && K[it1->first] > K[it2->first]) || (K[it1->second] > K[it2->second] && K[it1->first] < K[it2->first])) {
                            swap(it1->second,it2->second);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->first].size() == E[it2->second].size()) {
                        if ((K[it1->second] < K[it2->first] && K[it1->first] > K[it2->second]) || (K[it1->second] > K[it2->first] && K[it1->first] < K[it2->second])) {
                            swap(it1->second,it2->first);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->second].size() == E[it2->first].size()) {
                        if ((K[it1->first] < K[it2->second] && K[it1->second] > K[it2->first]) || (K[it1->first] > K[it2->second] && K[it1->second] < K[it2->first])) {
                            swap(it1->first,it2->second);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->second].size() == E[it2->second].size()) {
                        if ((K[it1->first] < K[it2->first] && K[it1->second] > K[it2->second]) || (K[it1->first] > K[it2->first] && K[it1->second] < K[it2->second])) {
                            swap(it1->first,it2->first);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                }
            }
        }
        double sum_kiki = 0;
        double sum_ki = 0;
        double sum_ki_2 = 0;
        for (auto e : edge_copy) {
            sum_kiki += K[e.first] * K[e.second];
            sum_ki += K[e.first] + K[e.second];
            sum_ki_2 += K[e.first] * K[e.first] + K[e.second] * K[e.second];
        }
        double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        wcout << L"Second order alpha: " << alpha << endl;
    }

    {
        // Link Rewiring (second order MAX disassortative)
        vector<double> K;
        K.resize(E.size());
        for (size_t i = 0 ; i < E.size() ; ++i) {
            K[i] = 0;
            for (auto u : E[i]) {
                K[i] += E[id[u.first]].size();
            }
            K[i] /= static_cast<double>(E[i].size());
        }
        vector<pair<size_t,size_t>> edge_copy;
        for (auto e : edge) {
            edge_copy.push_back(e);
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto it1 = edge_copy.begin() ; it1 != edge_copy.end() ; ++it1) {
                auto it2 = it1;
                for (++it2 ; it2 != edge_copy.end() ; ++it2) {
                    if (E[it1->first].size() == E[it2->first].size()) {
                        if ((K[it1->second] > K[it2->second] && K[it1->first] > K[it2->first]) || (K[it1->second] < K[it2->second] && K[it1->first] < K[it2->first])) {
                            swap(it1->second,it2->second);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->first].size() == E[it2->second].size()) {
                        if ((K[it1->second] > K[it2->first] && K[it1->first] > K[it2->second]) || (K[it1->second] < K[it2->first] && K[it1->first] < K[it2->second])) {
                            swap(it1->second,it2->first);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->second].size() == E[it2->first].size()) {
                        if ((K[it1->first] > K[it2->second] && K[it1->second] > K[it2->first]) || (K[it1->first] < K[it2->second] && K[it1->second] < K[it2->first])) {
                            swap(it1->first,it2->second);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->second].size() == E[it2->second].size()) {
                        if ((K[it1->first] > K[it2->first] && K[it1->second] > K[it2->second]) || (K[it1->first] < K[it2->first] && K[it1->second] < K[it2->second])) {
                            swap(it1->first,it2->first);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                }
            }
        }
        double sum_kiki = 0;
        double sum_ki = 0;
        double sum_ki_2 = 0;
        for (auto e : edge_copy) {
            sum_kiki += K[e.first] * K[e.second];
            sum_ki += K[e.first] + K[e.second];
            sum_ki_2 += K[e.first] * K[e.first] + K[e.second] * K[e.second];
        }
        double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        wcout << L"Second order alpha: " << alpha << endl;
    }

    {
        // Link Rewiring (second order MAX random)
        const size_t N = 30;
        default_random_engine generator(time(NULL));
        uniform_int_distribution<int> distribution(0,1);
        uniform_int_distribution<int> dis_2(0,edge.size() - 1);
        vector<double> K;
        K.resize(E.size());
        for (size_t i = 0 ; i < E.size() ; ++i) {
            K[i] = 0;
            for (auto u : E[i]) {
                K[i] += E[id[u.first]].size();
            }
            K[i] /= static_cast<double>(E[i].size());
        }
        double alpha_cum = 0;
        for (size_t rr = 0 ; rr < N ; ++rr) {
            vector<pair<size_t,size_t>> edge_copy;
            for (auto e : edge) {
                edge_copy.push_back(e);
            }

            size_t trial_constant = 10000 * total_number_link;
            for (size_t run = 0 ; run < trial_constant ; ++run) {
                auto it1 = &edge_copy[dis_2(generator)];
                auto it2 = &edge_copy[dis_2(generator)];
                if (E[it1->first].size() == E[it2->first].size() || E[it1->first].size() == E[it2->second].size() || E[it1->second].size() == E[it2->first].size() || E[it1->second].size() == E[it2->second].size()) {
                    size_t A[2];
                    size_t B[2];
                    A[0] = it1->first;
                    A[1] = it1->second;
                    B[0] = it2->first;
                    B[1] = it2->second;

                    size_t a0 = distribution(generator);
                    size_t b0 = distribution(generator);
                    while (E[A[a0]].size() != E[B[b0]].size()) {
                        a0 = distribution(generator);
                        b0 = distribution(generator);
                    }

                    size_t rewire = distribution(generator);
                    if (rewire) {
                        it1->first = A[a0];
                        it1->second = B[b0 ^ 1];
                        it2->first = B[b0];
                        it2->second = A[a0 ^ 1];
                    }
                }
            }

            double sum_kiki = 0;
            double sum_ki = 0;
            double sum_ki_2 = 0;
            for (auto e : edge_copy) {
                sum_kiki += K[e.first] * K[e.second];
                sum_ki += K[e.first] + K[e.second];
                sum_ki_2 += K[e.first] * K[e.first] + K[e.second] * K[e.second];
            }
            double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
            alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
            wcout << L"R alpha: " << alpha << endl;
            alpha_cum += alpha;
        }
        wcout << L"Avg alpha: " << alpha_cum / static_cast<double>(N) << endl;
    }

    {
        // Second order max assortative coef
        vector<double> K;
        K.resize(E.size());
        for (size_t i = 0 ; i < E.size() ; ++i) {
            K[i] = 0;
            for (auto u : E[i]) {
                if (K[i] < E[id[u.first]].size())
                    K[i] = E[id[u.first]].size();
            }
        }
        double sum_kiki = 0;
        double sum_ki = 0;
        double sum_ki_2 = 0;
        for (auto e : edge) {
            sum_kiki += K[e.first] * K[e.second];
            sum_ki += K[e.first] + K[e.second];
            sum_ki_2 += K[e.first] * K[e.first] + K[e.second] * K[e.second];
        }
        double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        wcout << L"Second order alpha: " << alpha << endl;
    }


    {
        // Link Rewiring (second order MAX assortative)
        vector<double> K;
        K.resize(E.size());
        for (size_t i = 0 ; i < E.size() ; ++i) {
            K[i] = 0;
            for (auto u : E[i]) {
                if (K[i] < E[id[u.first]].size())
                    K[i] = E[id[u.first]].size();
            }
        }
        vector<pair<size_t,size_t>> edge_copy;
        for (auto e : edge) {
            edge_copy.push_back(e);
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto it1 = edge_copy.begin() ; it1 != edge_copy.end() ; ++it1) {
                auto it2 = it1;
                for (++it2 ; it2 != edge_copy.end() ; ++it2) {
                    if (E[it1->first].size() == E[it2->first].size()) {
                        if ((K[it1->second] < K[it2->second] && K[it1->first] > K[it2->first]) || (K[it1->second] > K[it2->second] && K[it1->first] < K[it2->first])) {
                            swap(it1->second,it2->second);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->first].size() == E[it2->second].size()) {
                        if ((K[it1->second] < K[it2->first] && K[it1->first] > K[it2->second]) || (K[it1->second] > K[it2->first] && K[it1->first] < K[it2->second])) {
                            swap(it1->second,it2->first);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->second].size() == E[it2->first].size()) {
                        if ((K[it1->first] < K[it2->second] && K[it1->second] > K[it2->first]) || (K[it1->first] > K[it2->second] && K[it1->second] < K[it2->first])) {
                            swap(it1->first,it2->second);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->second].size() == E[it2->second].size()) {
                        if ((K[it1->first] < K[it2->first] && K[it1->second] > K[it2->second]) || (K[it1->first] > K[it2->first] && K[it1->second] < K[it2->second])) {
                            swap(it1->first,it2->first);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                }
            }
        }
        double sum_kiki = 0;
        double sum_ki = 0;
        double sum_ki_2 = 0;
        for (auto e : edge_copy) {
            sum_kiki += K[e.first] * K[e.second];
            sum_ki += K[e.first] + K[e.second];
            sum_ki_2 += K[e.first] * K[e.first] + K[e.second] * K[e.second];
        }
        double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        wcout << L"Second order alpha: " << alpha << endl;
    }

    {
        // Link Rewiring (second order MAX disassortative)
        vector<double> K;
        K.resize(E.size());
        for (size_t i = 0 ; i < E.size() ; ++i) {
            K[i] = 0;
            for (auto u : E[i]) {
                if (K[i] < E[id[u.first]].size())
                    K[i] = E[id[u.first]].size();
            }
        }
        vector<pair<size_t,size_t>> edge_copy;
        for (auto e : edge) {
            edge_copy.push_back(e);
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto it1 = edge_copy.begin() ; it1 != edge_copy.end() ; ++it1) {
                auto it2 = it1;
                for (++it2 ; it2 != edge_copy.end() ; ++it2) {
                    if (E[it1->first].size() == E[it2->first].size()) {
                        if ((K[it1->second] > K[it2->second] && K[it1->first] > K[it2->first]) || (K[it1->second] < K[it2->second] && K[it1->first] < K[it2->first])) {
                            swap(it1->second,it2->second);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->first].size() == E[it2->second].size()) {
                        if ((K[it1->second] > K[it2->first] && K[it1->first] > K[it2->second]) || (K[it1->second] < K[it2->first] && K[it1->first] < K[it2->second])) {
                            swap(it1->second,it2->first);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->second].size() == E[it2->first].size()) {
                        if ((K[it1->first] > K[it2->second] && K[it1->second] > K[it2->first]) || (K[it1->first] < K[it2->second] && K[it1->second] < K[it2->first])) {
                            swap(it1->first,it2->second);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                    if (E[it1->second].size() == E[it2->second].size()) {
                        if ((K[it1->first] > K[it2->first] && K[it1->second] > K[it2->second]) || (K[it1->first] < K[it2->first] && K[it1->second] < K[it2->second])) {
                            swap(it1->first,it2->first);
                            --it2;
                            changed = true;
                            continue;
                        }
                    }
                }
            }
        }
        double sum_kiki = 0;
        double sum_ki = 0;
        double sum_ki_2 = 0;
        for (auto e : edge_copy) {
            sum_kiki += K[e.first] * K[e.second];
            sum_ki += K[e.first] + K[e.second];
            sum_ki_2 += K[e.first] * K[e.first] + K[e.second] * K[e.second];
        }
        double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
        wcout << L"Second order alpha: " << alpha << endl;
    }

    {
        // Link Rewiring (second order MAX random)
        const size_t N = 30;
        default_random_engine generator(time(NULL));
        uniform_int_distribution<int> distribution(0,1);
        uniform_int_distribution<int> dis_2(0,edge.size() - 1);
        vector<double> K;
        K.resize(E.size());
        for (size_t i = 0 ; i < E.size() ; ++i) {
            K[i] = 0;
            for (auto u : E[i]) {
                if (K[i] < E[id[u.first]].size())
                    K[i] = E[id[u.first]].size();
            }
        }
        double alpha_cum = 0;
        for (size_t rr = 0 ; rr < N ; ++rr) {
            vector<pair<size_t,size_t>> edge_copy;
            for (auto e : edge) {
                edge_copy.push_back(e);
            }

            size_t trial_constant = 10000 * total_number_link;
            for (size_t run = 0 ; run < trial_constant ; ++run) {
                auto it1 = &edge_copy[dis_2(generator)];
                auto it2 = &edge_copy[dis_2(generator)];
                if (E[it1->first].size() == E[it2->first].size() || E[it1->first].size() == E[it2->second].size() || E[it1->second].size() == E[it2->first].size() || E[it1->second].size() == E[it2->second].size()) {
                    size_t A[2];
                    size_t B[2];
                    A[0] = it1->first;
                    A[1] = it1->second;
                    B[0] = it2->first;
                    B[1] = it2->second;

                    size_t a0 = distribution(generator);
                    size_t b0 = distribution(generator);
                    while (E[A[a0]].size() != E[B[b0]].size()) {
                        a0 = distribution(generator);
                        b0 = distribution(generator);
                    }

                    size_t rewire = distribution(generator);
                    if (rewire) {
                        it1->first = A[a0];
                        it1->second = B[b0 ^ 1];
                        it2->first = B[b0];
                        it2->second = A[a0 ^ 1];
                    }
                }
            }

            double sum_kiki = 0;
            double sum_ki = 0;
            double sum_ki_2 = 0;
            for (auto e : edge_copy) {
                sum_kiki += K[e.first] * K[e.second];
                sum_ki += K[e.first] + K[e.second];
                sum_ki_2 += K[e.first] * K[e.first] + K[e.second] * K[e.second];
            }
            double alpha = sum_kiki / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
            alpha /= sum_ki_2 / 2.0L / static_cast<double>(total_number_link) - (sum_ki / static_cast<double>(total_number_link) / 2.0L) * (sum_ki / static_cast<double>(total_number_link) / 2.0L);
            wcout << L"R alpha: " << alpha << endl;
            alpha_cum += alpha;
        }
        wcout << L"Avg alpha: " << alpha_cum / static_cast<double>(N) << endl;
    }
*/

/*
    {
        // Calculate betweenness value
        atomic<double> C_B[100000];
        for (size_t i = 0 ; i < 100000 ; ++i) {
            C_B[i] = 0;
        }
        auto f = [&](size_t start, size_t end) {
            for (size_t i = start ; i < end ; ++i) {
                //wcout << i << endl;
                queue<deque<size_t>> bfs;
                vector<vector<deque<size_t>>> dest_path;
                
                dest_path.resize(E.size());

                deque<size_t> temp;
                temp.push_back(i);
                bfs.push(temp);

                while (!bfs.empty()) {
                    temp = bfs.front();
                    bfs.pop();

                    size_t v = temp.back();

                    if (dest_path[v].empty()) {
                        dest_path[v].push_back(temp);

                        for (auto u : E[v]) {
                            temp.push_back(id[u.first]);
                            bfs.push(temp);
                            temp.pop_back();
                        }
                    } else if (dest_path[v][0].size() >= temp.size()) {
                        dest_path[v].push_back(temp);

                        for (auto u : E[v]) {
                            temp.push_back(id[u.first]);
                            bfs.push(temp);
                            temp.pop_back();
                        }
                    }
                }

                for (size_t j = 0 ; j < dest_path.size() ; ++j) {
                    double n = dest_path[j].size();
                    for (auto l : dest_path[j]) {
                        for (size_t k = 1 ; k < l.size() - 1 ; ++k) {
                            double temp = C_B[l[k]];
                            while (!C_B[l[k]].compare_exchange_strong(temp,temp + static_cast<double>(1.0L / n))) {}
                            //C_B[l[k]] += static_cast<double>(1.0L / n);
                        }
                    }
                }
            }
        };
        vector<thread> thread_pool;
        for (size_t i = 0 ; i < 8 ; ++i) {
            thread_pool.push_back(std::thread(f,i * E.size() / 8, (i + 1) * E.size() / 8));
        }
        for (size_t i = 0 ; i < 8 ; ++i) {
            thread_pool[i].join();
        }

        for (size_t i = 0 ; i < E.size() ; ++i) {
            wcout << r_id[i] << L"\t" << C_B[i] / 2.0L << endl;
        }
    }
*/
/*
    {
        // Generate BOSAM matrix
        vector<vector<size_t>> data;
        data.resize(E.size());
        for (size_t i = 0 ; i < E.size() ; ++i) {
            data[i].resize(4);
            data[i][0] = E[i].size();
            for (auto v : E[i]) {
                if (data[i][1] < E[id[v.first]].size()) {
                    data[i][1] = E[id[v.first]].size();
                }
                if (id[v.first] > data[i][2]) {
                    data[i][2] = id[v.first];
                }
            }
            data[i][3] = i;
        }
        sort(data.begin(),data.end());
        vector<vector<size_t>> data_copy = data;
        for (size_t i = 0 ; i < data.size() ; ++i) {
            for (size_t j = 0 ; j < 4 ; ++j) {
                data[i][j] = data_copy[data.size() - 1 - i][j];
            }
        }
        vector<size_t> reverse;
        reverse.resize(E.size());
        count = E.size();
        for (size_t i = 0 ; i < E.size() ; ++i) {
            reverse[data[i][3]] = i;
            if (data[i][0] <= 0) {
                count = i;
                break;
            }
        }

        vector<vector<bool>> output;
        for (size_t i = 0 ; i < count ; ++i) {
            vector<bool> temp;
            temp.resize(count);
            for (auto e : E[data[i][3]]) {
                temp[reverse[id[e.first]]] = true;
            }
            output.push_back(temp);
        }

        for (auto a : output) {
            for (auto b : a) {
                wcout << b << L" ";
            }
            wcout << endl;
        }
    }
*/
/*
    // Degree distribution
    count = 1;
    for (size_t i = 1 ; i <= Esize.size() ; ++i) {
        if (i == Esize.size() || Esize[i].first != Esize[i - 1].first) {
            wcout << L"(" << Esize[i-1].first << L"," << static_cast<double>(count) / static_cast<double>(Esize.size()) << L")" << endl;
            count = 1;
        } else {
            ++count;
        }
    }
*/
    return 0;
}