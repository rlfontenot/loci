std::vector<char>   extract_general_face(const Entity* lower, int lower_size,
                                         const Entity* upper, int upper_size,
                                         const Entity* boundary_map, int boundary_map_size,
                                         const const_multiMap& face2node,
                                         const const_multiMap& face2edge,
                                         const const_MapVec<2>& edge2node,
                                         const std::vector<char>& cellPlan,
                                         Entity ff,
                                         const Map& node_remap
                                         );
std::vector<char>  extract_prism_face(const  std::vector<char>& cellPlan, int dd);
std::vector<char>  merge_tri_face_p(const  std::vector<char>& cellPlan1, int dd1, char orientCode1);
