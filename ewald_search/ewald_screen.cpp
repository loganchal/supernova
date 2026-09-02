#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;
using Vec = std::array<int8_t,8>;
constexpr int D = 8;
constexpr int NPTS = 6561;
constexpr int WORDS = (NPTS + 63) / 64;

struct VecHash {
    std::size_t operator()(Vec const& a) const noexcept {
        std::uint64_t h = 1469598103934665603ULL;
        for (auto x : a) { h ^= static_cast<std::uint8_t>(x); h *= 1099511628211ULL; }
        return static_cast<std::size_t>(h);
    }
};

static std::vector<Vec> points;
static std::array<std::array<int,256>,9> subset_index{};
static std::array<std::vector<int>,9> subsets;
static std::unordered_map<Vec, std::vector<std::uint64_t>, VecHash> valid_cache;
static long long cache_hits=0, cache_misses=0;
static int global_min_coord=0, global_max_coord=0;

static Vec canonical_q(Vec q) {
    for (int i=0;i<D;i++) if (q[i]) {
        if (q[i] < 0) for (auto &z:q) z = -z;
        break;
    }
    return q;
}

static const std::vector<std::uint64_t>& validity_mask(Vec q) {
    q = canonical_q(q);
    auto it=valid_cache.find(q);
    if (it!=valid_cache.end()) { ++cache_hits; return it->second; }
    ++cache_misses;
    std::vector<std::uint64_t> m(WORDS,0);
    for (int i=0;i<NPTS;i++) {
        int s=0;
        for (int j=0;j<D;j++) s += int(q[j])*int(points[i][j]);
        if (-1 <= s && s <= 1) m[i>>6] |= (1ULL<<(i&63));
    }
    return valid_cache.emplace(q,std::move(m)).first->second;
}

static long long det_matrix(std::array<std::array<long long,D>,D> a) {
    long long sign=1, prev=1;
    for (int k=0;k<D-1;k++) {
        int p=k;
        while (p<D && a[p][k]==0) ++p;
        if (p==D) return 0;
        if (p!=k) { std::swap(a[p],a[k]); sign=-sign; }
        long long pivot=a[k][k];
        for (int i=k+1;i<D;i++) for (int j=k+1;j<D;j++) {
            __int128 num=(__int128)a[i][j]*pivot-(__int128)a[i][k]*a[k][j];
            if (k>0) num/=prev;
            a[i][j]=(long long)num;
        }
        prev=pivot;
    }
    return sign*a[D-1][D-1];
}

static long long det_basis(const std::vector<Vec>& s, const std::array<int,D>& idx) {
    std::array<std::array<long long,D>,D> a{};
    for (int i=0;i<D;i++) for(int j=0;j<D;j++) a[i][j]=s[idx[i]][j];
    return det_matrix(a);
}

struct F2Basis { std::array<uint8_t,D> piv{}; int rank=0; };
static bool f2_insert(F2Basis& b, uint8_t x) {
    for (int bit=D-1;bit>=0;--bit) if ((x>>bit)&1) {
        if (b.piv[bit]) x ^= b.piv[bit];
        else {
            b.piv[bit]=x; ++b.rank;
            for (int j=0;j<D;j++) if (j!=bit && ((b.piv[j]>>bit)&1)) b.piv[j]^=x;
            return true;
        }
    }
    return false;
}

static int rank_mod2(const std::vector<Vec>& s) {
    F2Basis b;
    for (auto const& v:s) {
        uint8_t x=0; for(int j=0;j<D;j++) if(v[j]&1) x|=(1u<<j);
        f2_insert(b,x);
    }
    return b.rank;
}

static std::vector<long long> wedge_extend(const std::vector<long long>& w, int k, const Vec& v) {
    std::vector<long long> out(subsets[k+1].size());
    for (std::size_t z=0;z<subsets[k+1].size();++z) {
        int mask=subsets[k+1][z], t=0;
        long long sum=0;
        for (int j=0;j<D;j++) if ((mask>>j)&1) {
            int sub=mask^(1<<j);
            long long term=(long long)v[j]*w[subset_index[k][sub]];
            if ((k+t)&1) sum-=term; else sum+=term;
            ++t;
        }
        out[z]=sum;
    }
    return out;
}

static long long wedge_gcd(const std::vector<long long>& w) {
    long long g=0;
    for (auto x:w) g=std::gcd(g,std::llabs(x));
    return g;
}

struct SearchResult {
    bool found=false;
    bool capped=false;
    uint64_t nodes=0;
    long long min_abs_det=0;
    std::array<int,D> basis{};
};

static bool dfs_basis(const std::vector<Vec>& s, int start, int k,
                      const std::vector<long long>& w, F2Basis f2,
                      std::array<int,D>& chosen, SearchResult& r, uint64_t cap) {
    if (r.nodes++ >= cap) { r.capped=true; return false; }
    if ((int)s.size()-start < D-k) return false;
    if (k==D) {
        long long d=std::llabs(w[0]);
        if (!r.min_abs_det || d<r.min_abs_det) r.min_abs_det=d;
        if (d==1) { r.found=true; r.basis=chosen; return true; }
        return false;
    }
    for (int i=start;i<=(int)s.size()-(D-k);++i) {
        F2Basis nf=f2;
        uint8_t pm=0; for(int j=0;j<D;j++) if(s[i][j]&1) pm|=(1u<<j);
        if (!f2_insert(nf,pm)) continue;
        auto nw=wedge_extend(w,k,s[i]);
        if (wedge_gcd(nw)!=1) continue;
        chosen[k]=i;
        if (dfs_basis(s,i+1,k+1,nw,nf,chosen,r,cap)) return true;
        if (r.capped) return false;
    }
    return false;
}

static SearchResult find_unimodular_basis(std::vector<Vec>& s, uint64_t cap, uint64_t seed) {
    SearchResult r;
    if (s.size()<D || rank_mod2(s)<D) return r;
    auto support=[](Vec const& v){int z=0;for(auto x:v)z+=x!=0;return z;};
    std::stable_sort(s.begin(),s.end(),[&](Vec const&a,Vec const&b){
        int sa=support(a), sb=support(b); if(sa!=sb)return sa<sb; return a<b;
    });

    std::vector<int> order(s.size()); std::iota(order.begin(),order.end(),0);
    std::mt19937_64 rng(seed);
    int trials=std::min<int>(64, 4+(int)s.size());
    for(int tr=0;tr<trials;tr++) {
        if(tr>0) std::shuffle(order.begin(),order.end(),rng);
        F2Basis f; std::array<int,D> idx{}; int k=0;
        for(int i:order) {
            F2Basis nf=f; uint8_t pm=0;
            for(int j=0;j<D;j++)if(s[i][j]&1)pm|=1u<<j;
            if(f2_insert(nf,pm)){f=nf;idx[k++]=i;if(k==D)break;}
        }
        if(k==D){
            long long d=std::llabs(det_basis(s,idx));
            if(!r.min_abs_det||d<r.min_abs_det)r.min_abs_det=d;
            if(d==1){r.found=true;r.basis=idx;return r;}
        }
    }

    std::array<int,D> chosen{};
    std::vector<long long> w0{1}; F2Basis f0;
    dfs_basis(s,0,0,w0,f0,chosen,r,cap);
    return r;
}

static std::vector<Vec> decode_line(const std::string& line, int base) {
    cpp_int n(line);
    std::vector<int> a;
    while(n>0) {
        cpp_int q=n/base; cpp_int r=n-q*base;
        a.push_back(r.convert_to<int>()); n=q;
    }
    if(a.size()<2) throw std::runtime_error("short encoded record");
    int dim=a[0], shift=a[1];
    if(dim!=D || (a.size()-2)%D) throw std::runtime_error("bad encoded record dimensions");
    std::vector<Vec> vs;
    for(std::size_t i=2;i<a.size();i+=D){
        Vec v{};
        for(int j=0;j<D;j++){
            int x=a[i+j]-shift;
            if(x<-127||x>127)throw std::runtime_error("coordinate outside int8 range");
            v[j]=(int8_t)x;
            global_min_coord=std::min(global_min_coord,x);
            global_max_coord=std::max(global_max_coord,x);
        }
        vs.push_back(v);
    }
    return vs;
}

static bool has_standard_units(const std::vector<Vec>& q) {
    for(int j=0;j<D;j++){
        bool ok=false;
        for(auto const&v:q){
            bool same=true;
            for(int k=0;k<D;k++)if(v[k]!=(k==j)){same=false;break;}
            if(same){ok=true;break;}
        }
        if(!ok)return false;
    }
    return true;
}

static std::string vec_json(Vec const& v) {
    std::ostringstream o; o<<"[";
    for(int j=0;j<D;j++){if(j)o<<",";o<<int(v[j]);}
    o<<"]"; return o.str();
}
static std::string vecs_json(std::vector<Vec> const& vs) {
    std::ostringstream o; o<<"[";
    for(std::size_t i=0;i<vs.size();i++){if(i)o<<",";o<<vec_json(vs[i]);}
    o<<"]"; return o.str();
}

int main(int argc,char**argv){
    if(argc<4){std::cerr<<"usage: screen DATA_DIR BLOCK_LO BLOCK_HI [CAP]\n";return 2;}
    std::string dir=argv[1]; int blo=std::stoi(argv[2]), bhi=std::stoi(argv[3]);
    uint64_t cap=argc>=5?std::stoull(argv[4]):200000;
    for(int k=0;k<=D;k++){
        subset_index[k].fill(-1);
        for(int m=0;m<256;m++)if(std::popcount((unsigned)m)==k){
            subset_index[k][m]=subsets[k].size(); subsets[k].push_back(m);
        }
    }
    points.reserve(NPTS);
    for(int code=0;code<NPTS;code++){
        int z=code; Vec v{};
        for(int j=0;j<D;j++){int dig=z%3;z/=3;v[j]=(int8_t)(dig-1);}
        points.push_back(v);
    }

    std::uint64_t total=0,basis_found=0,unresolved=0,rankdef=0,bad_units=0;
    std::uint64_t min_ecount=UINT64_MAX,max_ecount=0;
    auto t0=std::chrono::steady_clock::now();
    std::cout<<"{\"type\":\"header\",\"blocks\":["<<blo<<","<<bhi<<"],\"cap\":"<<cap<<"}\n";
    for(int b=blo;b<=bhi;b++){
        std::ifstream f(dir+"/block"+std::to_string(b));
        if(!f){std::cerr<<"missing block "<<b<<"\n";continue;}
        std::string line; std::getline(f,line); int base=std::stoi(line),rec=0;
        while(std::getline(f,line)){
            if(line.empty())continue; ++rec; ++total;
            auto q=decode_line(line,base);
            int id=b*7498+rec;
            if(!has_standard_units(q)){
                ++bad_units;
                std::cout<<"{\"type\":\"bad_units\",\"id\":"<<id<<",\"normals\":"<<vecs_json(q)<<"}\n";
                continue;
            }
            std::vector<uint64_t> em(WORDS,~0ULL);
            if(NPTS%64)em.back()&=((1ULL<<(NPTS%64))-1);
            for(auto const&v:q){
                auto const&m=validity_mask(v);
                for(int w=0;w<WORDS;w++)em[w]&=m[w];
            }
            std::uint64_t ec=0; for(auto x:em)ec+=std::popcount(x);
            min_ecount=std::min(min_ecount,ec); max_ecount=std::max(max_ecount,ec);
            std::vector<Vec> ehalf;
            for(int i=0;i<NPTS;i++)if((em[i>>6]>>(i&63))&1ULL){
                auto v=points[i]; bool zero=true,pos=false;
                for(int j=0;j<D;j++)if(v[j]){zero=false;pos=v[j]>0;break;}
                if(!zero&&pos)ehalf.push_back(v);
            }
            int r2=rank_mod2(ehalf);
            if(r2<D){
                ++rankdef;
                std::cout<<"{\"type\":\"rank2_counterexample\",\"id\":"<<id<<",\"block\":"<<b<<",\"record\":"<<rec<<",\"nvertices\":"<<q.size()<<",\"ecount\":"<<ec<<",\"rank2\":"<<r2<<",\"encoded\":\""<<line<<"\",\"normals\":"<<vecs_json(q)<<",\"ehalf\":"<<vecs_json(ehalf)<<"}\n";
                continue;
            }
            auto sr=find_unimodular_basis(ehalf,cap,0x9E3779B97F4A7C15ULL^(uint64_t)id);
            if(sr.found)++basis_found;
            else{
                ++unresolved;
                std::cout<<"{\"type\":\"unresolved\",\"id\":"<<id<<",\"block\":"<<b<<",\"record\":"<<rec<<",\"nvertices\":"<<q.size()<<",\"ecount\":"<<ec<<",\"rank2\":"<<r2<<",\"nodes\":"<<sr.nodes<<",\"capped\":"<<(sr.capped?"true":"false")<<",\"min_abs_det_seen\":"<<sr.min_abs_det<<",\"encoded\":\""<<line<<"\",\"normals\":"<<vecs_json(q)<<",\"ehalf\":"<<vecs_json(ehalf)<<"}\n";
            }
            if(total%10000==0){
                double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
                std::cerr<<"processed="<<total<<" found="<<basis_found<<" unresolved="<<unresolved<<" rankdef="<<rankdef<<" cache="<<valid_cache.size()<<" sec="<<sec<<"\n";
            }
        }
    }
    double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::cout<<"{\"type\":\"summary\",\"total\":"<<total<<",\"basis_found\":"<<basis_found<<",\"unresolved\":"<<unresolved<<",\"rank2_counterexamples\":"<<rankdef<<",\"bad_units\":"<<bad_units<<",\"min_ecount\":"<<(min_ecount==UINT64_MAX?0:min_ecount)<<",\"max_ecount\":"<<max_ecount<<",\"cache_size\":"<<valid_cache.size()<<",\"cache_hits\":"<<cache_hits<<",\"cache_misses\":"<<cache_misses<<",\"coord_min\":"<<global_min_coord<<",\"coord_max\":"<<global_max_coord<<",\"seconds\":"<<std::setprecision(8)<<sec<<"}\n";
}
