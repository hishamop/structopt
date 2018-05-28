#ifndef MATERIALMANAGER_H
#define MATERIALMANAGER_H
#include"material/cmaterial.h"
#include<memory>

using material_ptr= std::shared_ptr<CMaterial>;
class MaterialManager
{
public:
    MaterialManager();
    static material_ptr get_material();

private:
    static material_ptr m_material;
};

#endif // MATERIALMANAGER_H
