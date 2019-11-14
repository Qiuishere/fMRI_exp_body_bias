function List = listdir(path_to_dir)

List = dir(path_to_dir);
List = {List.name}';

List(strncmp('.',List,1)) = [];

end