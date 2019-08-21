// python stuff

static inline void initPython() {
    static bool must_init=true;
    if (must_init) {
         Py_Initialize(); // init boost & numpy boost
         np::initialize();
         must_init=false;
    }
}

template <class T>
void*clone_data(vector<T>&v) { //  vector data clone
    auto sz=v.size()*sizeof(T);
    return memcpy(malloc(sz), v.data(), sz);
}

// create a cloned data numpy array
template <class T>
static np::ndarray vector2numpy(vector<T>v) {
    initPython();
    return np::from_data(clone_data(v),     // data -> clone
            np::dtype::get_builtin<T>(),  // dtype -> double
            p::make_tuple(v.size()),    // shape -> size
            p::make_tuple(sizeof(T)), p::object()); // stride
}

template <class T>
static np::ndarray vector2Coords(vector<T>v) {
    initPython();
    return np::from_data(clone_data(v),     // data -> clone
            np::dtype::get_builtin<T>(),  // dtype -> double
            p::make_tuple(v.size()/3, 3),    // shape
            p::make_tuple(sizeof(T)*3, sizeof(T)), p::object()); // stride
}

template <class T>
static np::ndarray vector2Textures(vector<T>v) {
    initPython();
    return np::from_data(clone_data(v),     // data -> clone
            np::dtype::get_builtin<T>(),  // dtype -> double
            p::make_tuple(v.size()/2, 2),    // shape
            p::make_tuple(sizeof(T)*2, sizeof(T)), p::object()); // stride
}

template <class T>
vector<T>  normalizeCoords(vector<T>coords) { // 0..1
    auto   mm = std::minmax_element(coords.begin(), coords.end());
    double diff=abs(*mm.second - *mm.first);
    if(diff!=0)
        for (size_t i=0; i<coords.size(); i++) coords[i]/=diff;
    return coords;
}

template <class T>
p::list vector2list(vector<T>v) {
    p::list l;
    for (auto _v:v) l.append(_v);
    return l;
}

template <class T>
p::list vv2list(vector<T>v) { // vector of vector to list
    p::list l;
    for (auto _v:v)
        l.append( vector2list(_v) );
    return l;
}

vector<float>list2vector(p::list l) {
    vector<float>v;
    for (auto i=0; i<len(l); i++)
        v.push_back(p::extract<double>(l[i]));
    return v;
}

p::list mesh2list(Mesh mesh) {
    p::list l;
    l.append( vector2numpy    ( std::get<0>(mesh) ) ); //  indexes, coords, textures, normals
    l.append( vector2Coords   ( std::get<1>(mesh) ) );
    l.append( vector2Textures ( std::get<2>(mesh) ) );
    l.append( vector2Coords   ( std::get<3>(mesh) ) );
    return l;
}
