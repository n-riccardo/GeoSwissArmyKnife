const global where_my_CPTS_are = @__DIR__

function getMyCPTsCollectionPath()
    return where_my_CPTS_are
end

function getMyCPTsList()
    return readdir(getMyCPTsCollectionPath())
end

function getMyCPTPath(cpt_name)
    # Get the path only if the CPT exists
    if cpt_name in getMyCPTsList()
        return getMyCPTsCollectionPath() * "/" * cpt_name
    else
        print("Your cpt does not exist in path!")
        return nothing
    end
end

function lookCPTcontent(cpt_name)
    # Get the path only if the CPT exists
    if cpt_name in getMyCPTsList()
        return readlines(getMyCPTPath(cpt_name))
    else
        print("Your cpt does not exist in path!")
        return nothing
    end
end