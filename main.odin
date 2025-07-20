package main
import "core:fmt"
import "core:testing"
import "core:math"
import "core:thread"
import "core:strings"
import "vendor:raylib"

SCALE :: 100
WIDTH :: 400*2
HEIGHT :: 250*2
log :: proc(a: ..any){}//fmt.print
logln :: proc(a: ..any){}//fmt.println


// get coords to draw on screen based on center 
get_draw_coords :: proc(p :Vec2) -> Vec2 {
    return Vec2{WIDTH / 2 + p.x* SCALE, HEIGHT/2+p.y*SCALE}
}

Vec2 :: raylib.Vector2
sin::math.sin
cos::math.cos

normalize2 :: proc(v: Vec2) -> Vec2 {
    mag := math.sqrt(v.x*v.x + v.y*v.y);
    if mag == 0.0 { return Vec2{0, 0}; }
    return Vec2{v.x / mag, v.y / mag};
}

normalize3 :: proc(v: Vec3) -> Vec3 {
    mag := math.sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if mag == 0.0 { return Vec3{0, 0, 0}; }
    return Vec3{v.x / mag, v.y / mag, v.z / mag};
}
rotate_point :: proc(p, origin: Vec2, angle_rad: f32) -> Vec2 {
    dx := p.x - origin.x;
    dy := p.y - origin.y;

    cos_a := cos(angle_rad);
    sin_a := sin(angle_rad);

    return Vec2{
        origin.x + dx * cos_a - dy * sin_a,
        origin.y + dx * sin_a + dy * cos_a,
    };
}

rotate_rectangle :: proc(points: [4]Vec2, origin: Vec2, angle_rad: f32) -> [4]Vec2 {
    result: [4]Vec2;
    for i in 0..<4 {
        result[i] = rotate_point(points[i], origin, angle_rad);
    }
    return result;
}
// METHOD 2: Simplified version (easier to understand)
ray_box_intersection_simple :: proc(ray: Ray, bbox: BoundingBox) -> bool {
    tmin : f32 = 0.0
    tmax : f32 = math.inf_f32(1)
    
    // Test intersection with each pair of parallel planes
    for axis in 0..<3 {
        ray_origin := axis == 0 ? ray.position.x : (axis == 1 ? ray.position.y : ray.position.z)
        ray_dir := axis == 0 ? ray.direction.x : (axis == 1 ? ray.direction.y : ray.direction.z)
        box_min := axis == 0 ? bbox.min.x : (axis == 1 ? bbox.min.y : bbox.min.z)
        box_max := axis == 0 ? bbox.max.x : (axis == 1 ? bbox.max.y : bbox.max.z)
        
        if abs(ray_dir) < 1e-8 {
            // Ray is parallel to this pair of planes
            if ray_origin < box_min || ray_origin > box_max {
                return false  // Ray is outside the slab
            }
        } else {
            // Calculate intersection distances
            t1 := (box_min - ray_origin) / ray_dir
            t2 := (box_max - ray_origin) / ray_dir
            
            // Make sure t1 is the near intersection
            if t1 > t2 {
                t1, t2 = t2, t1
            }
            
            // Update tmin and tmax
            tmin = max(tmin, t1)
            tmax = min(tmax, t2)
            
            // Early exit if no intersection possible
            if tmin > tmax {
                return false
            }
        }
    }
    
    return tmax >= 0  // Intersection exists and is in front of ray
}

Rect::raylib.Rectangle
Vec3::raylib.Vector3
Color::raylib.Color
Ray::raylib.Ray
Body :: struct {
    body: Rect,
    position: Vec2,
    rotation: f64,
    mass:   f64,
    color: Color,

    velocity: f64,
    angular_velocity:f64,

    acceleration:f64,
    angular_acceleration:f64,

    force: Vec2,
    torque: f64,

    inertia: f64,

    // coefficion of restetution
    cor: f64,
    // coefficion of friction
    cof: f64,

    mesh: []Triangle,
    bbox: BoundingBox,
}
BoundingBox :: struct {
    min: Vec3,  // Bottom-left-back corner
    max: Vec3,  // Top-right-front corner
}

dot :: proc(a, b: Vec3) -> f32 {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
cross :: proc(a, b: Vec3) -> Vec3 {
    return Vec3{
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x,
    };
}

Triangle :: [3]Vec3
// Initialize all your bodies with bounding boxes
init_bodies_with_bbox :: proc(bodies: []Body) -> []Body {
    result := make([]Body, len(bodies))
    
    for body, i in bodies {
        result[i] = body
        result[i].bbox = BoundingBox{}
    }
    
    return result
}

MollerTrumbore :: proc(v1, v2, v0 :Vec3, ray:Ray) -> (bool, f32) {
    e1 := v1 - v0;
    e2 := v2 - v0;

    p := cross (ray.direction, e2);
    det := dot(p, e1);
    epsilon : f32 = 1e-2;
    if (abs(det) <= epsilon) {
        return false, 0;
    }
    inv_det := 1.0/det;


    t_vec := ray.position - v0;

    u := (dot(t_vec,p)) * inv_det;
    if u < 0|| u > 1 {
        return false,0;
    }

    q := cross(t_vec,e1);

    v := dot(ray.direction, q)*inv_det;
    if v < 0|| v > 1 {
        return false,0;
    }
    if u + v > 1 {
        return false, 0;
    }

    t := dot(e2, q) * inv_det;
    if t > 0 {
        return true, t;
    }

    return false, 0;
}


draw_body :: proc(b: Body) {
    raylib.DrawRectanglePro(
        b.body, b.position, f32(b.rotation), b.color
    )
}

a : quaternion256

WORLD_UP :: Vec3{0,1,0};
Camera ::struct {
    pos: Vec3,
    c_forward: Vec3,
    focal_length: f32,
    c_right: Vec3,
    c_up: Vec3,
    s_center: Vec3, // c_forward * focal_length
    s_width: f32,
    s_height: f32,
    s_horizontal: Vec3, // c_right * s_width
    s_vertical: Vec3, // c_up * s_height

    c_angle_h: f32,
    c_angle_v: f32,
}
update_camera:: proc(c: ^Camera)-> Camera {
    s_ah := sin(c.c_angle_h)
    c_ah := cos(c.c_angle_h)
    s_av := sin(c.c_angle_v)
    c_av := cos(c.c_angle_v)
    d_ := c_av
    c.c_forward = normalize3(Vec3{
        -s_ah * c_av,
        s_av,
        c_ah * c_av,
    });

    c.s_center = c.c_forward * c.focal_length;
    c.c_right = normalize3(cross(WORLD_UP, c.c_forward));
    c.c_up = normalize3(cross(c.c_right,c.c_forward));

    c.s_horizontal = c.c_right * c.s_width
    c.s_vertical = c.c_up * c.s_height
    return c^
}
camera_bottom_left :: proc(c: ^Camera) -> Vec3 {
    return c.s_center - c.s_horizontal * 0.5 - c.s_vertical * 0.5;
}

camera_bottom_right :: proc(c: ^Camera) -> Vec3 {
    return c.s_center + c.s_horizontal * 0.5 - c.s_vertical * 0.5;
}

camera_top_left :: proc(c: ^Camera) -> Vec3 {
    return c.s_center - c.s_horizontal * 0.5 + c.s_vertical * 0.5;
}

camera_top_right :: proc(c: ^Camera) -> Vec3 {
    return c.s_center + c.s_horizontal * 0.5 + c.s_vertical * 0.5;
}


main :: proc() {
    logln("Hello, Oding")
    raylib.InitWindow(WIDTH, HEIGHT, "Engine 1")

    c:= Camera {
        Vec3{0,0,0},
        Vec3{0,0,1},
        10,
        Vec3{0,0,0},
        Vec3{0,0,0},
        Vec3{0,0,0},
        WIDTH/SCALE,
        HEIGHT/SCALE,
        Vec3{0,0,0},
        Vec3{0,0,0},
        0,
        0,
    };
    c =update_camera(&c);
    logln(c)
    bodies := init_bodies_with_bbox([]Body{
        Body{
            mesh= []Triangle{
                // CUBE 1 - Center cube (your original)
                // Front face (z = +0.05)
                {Vec3{-0.5, -0.5, 0.5}, Vec3{0.5, -0.5, 0.5}, Vec3{0.5, 0.5, 0.5}},
                {Vec3{-0.5, -0.5, 0.5}, Vec3{0.5, 0.5, 0.5}, Vec3{-0.5, 0.5, 0.5}},
                // Back face (z = -0.5)
                {Vec3{0.5, -0.5, -0.5}, Vec3{-0.5, -0.5, -0.5}, Vec3{-0.5, 0.5, -0.5}},
                {Vec3{0.5, -0.5, -0.5}, Vec3{-0.5, 0.5, -0.5}, Vec3{0.5, 0.5, -0.5}},
                // Left face (x = -0.5)
                {Vec3{-0.5, -0.5, -0.5}, Vec3{-0.5, -0.5, 0.5}, Vec3{-0.5, 0.5, 0.5}},
                {Vec3{-0.5, -0.5, -0.5}, Vec3{-0.5, 0.5, 0.5}, Vec3{-0.5, 0.5, -0.5}},
                // Right face (x = +0.5)
                {Vec3{0.5, -0.5, 0.5}, Vec3{0.5, -0.5, -0.5}, Vec3{0.5, 0.5, -0.5}},
                {Vec3{0.5, -0.5, 0.5}, Vec3{0.5, 0.5, -0.5}, Vec3{0.5, 0.5, 0.5}},
                // Top face (y = +0.5)
                {Vec3{-0.5, 0.5, 0.5}, Vec3{0.5, 0.5, 0.5}, Vec3{0.5, 0.5, -0.5}},
                {Vec3{-0.5, 0.5, 0.5}, Vec3{0.5, 0.5, -0.5}, Vec3{-0.5, 0.5, -0.5}},
                // Bottom face (y = -0.5)
                {Vec3{-0.5, -0.5, -0.5}, Vec3{0.5, -0.5, -0.5}, Vec3{0.5, -0.5, 0.5}},
                {Vec3{-0.5, -0.5, -0.5}, Vec3{0.5, -0.5, 0.5}, Vec3{-0.5, -0.5, 0.5}},
            },
            color=raylib.RED,
        },/*
        Body{
            mesh= []Triangle{

                // CUBE 2 - Left cube (offset by -1.0 in X)
                // Front face
                {Vec3{-1.5, -0.5, 0.5}, Vec3{-0.5, -0.5, 0.5}, Vec3{-0.5, 0.5, 0.5}},
                {Vec3{-1.5, -0.5, 0.5}, Vec3{-0.5, 0.5, 0.5}, Vec3{-1.5, 0.5, 0.5}},
                // Back face
                {Vec3{-0.5, -0.5, -0.5}, Vec3{-1.5, -0.5, -0.5}, Vec3{-1.5, 0.5, -0.5}},
                {Vec3{-0.5, -0.5, -0.5}, Vec3{-1.5, 0.5, -0.5}, Vec3{-0.5, 0.5, -0.5}},
                // Left face
                {Vec3{-1.5, -0.5, -0.5}, Vec3{-1.5, -0.5, 0.5}, Vec3{-1.5, 0.5, 0.5}},
                {Vec3{-1.5, -0.5, -0.5}, Vec3{-1.5, 0.5, 0.5}, Vec3{-1.5, 0.5, -0.5}},
                // Right face
                {Vec3{-0.5, -0.5, 0.5}, Vec3{-0.5, -0.5, -0.5}, Vec3{-0.5, 0.5, -0.5}},
                {Vec3{-0.5, -0.5, 0.5}, Vec3{-0.5, 0.5, -0.5}, Vec3{-0.5, 0.5, 0.5}},
                // Top face
                {Vec3{-1.5, 0.5, 0.5}, Vec3{-0.5, 0.5, 0.5}, Vec3{-0.5, 0.5, -0.5}},
                {Vec3{-1.5, 0.5, 0.5}, Vec3{-0.5, 0.5, -0.5}, Vec3{-1.5, 0.5, -0.5}},
                // Bottom face
                {Vec3{-1.5, -0.5, -0.5}, Vec3{-0.5, -0.5, -0.5}, Vec3{-0.5, -0.5, 0.5}},
                {Vec3{-1.5, -0.5, -0.5}, Vec3{-0.5, -0.5, 0.5}, Vec3{-1.5, -0.5, 0.5}},
            },color=raylib.BLUE,
        },
        Body{
            mesh= []Triangle{

                // CUBE 3 - Right cube (offset by +1.0 in X)
                // Front face
                {Vec3{0.5, -0.5, 0.5}, Vec3{1.5, -0.5, 0.5}, Vec3{1.5, 0.5, 0.5}},
                {Vec3{0.5, -0.5, 0.5}, Vec3{1.5, 0.5, 0.5}, Vec3{0.5, 0.5, 0.5}},
                // Back face
                {Vec3{1.5, -0.5, -0.5}, Vec3{0.5, -0.5, -0.5}, Vec3{0.5, 0.5, -0.5}},
                {Vec3{1.5, -0.5, -0.5}, Vec3{0.5, 0.5, -0.5}, Vec3{1.5, 0.5, -0.5}},
                // Left face
                {Vec3{0.5, -0.5, -0.5}, Vec3{0.5, -0.5, 0.5}, Vec3{0.5, 0.5, 0.5}},
                {Vec3{0.5, -0.5, -0.5}, Vec3{0.5, 0.5, 0.5}, Vec3{0.5, 0.5, -0.5}},
                // Right face
                {Vec3{1.5, -0.5, 0.5}, Vec3{1.5, -0.5, -0.5}, Vec3{1.5, 0.5, -0.5}},
                {Vec3{1.5, -0.5, 0.5}, Vec3{1.5, 0.5, -0.5}, Vec3{1.5, 0.5, 0.5}},
                // Top face
                {Vec3{0.5, 0.5, 0.5}, Vec3{1.5, 0.5, 0.5}, Vec3{1.5, 0.5, -0.5}},
                {Vec3{0.5, 0.5, 0.5}, Vec3{1.5, 0.5, -0.5}, Vec3{0.5, 0.5, -0.5}},
                // Bottom face
                {Vec3{0.5, -0.5, -0.5}, Vec3{1.5, -0.5, -0.5}, Vec3{1.5, -0.5, 0.5}},
                {Vec3{0.5, -0.5, -0.5}, Vec3{1.5, -0.5, 0.5}, Vec3{0.5, -0.5, 0.5}},
            },color=raylib.GREEN,
        },
        Body{
            mesh=[]Triangle{

                // CUBE 4 - Top cube (offset by +1.0 in Y)
                // Front face
                {Vec3{-05, 0.5, 05}, Vec3{05, 0.5, 05}, Vec3{05, 15, 05}},
                {Vec3{-05, 0.5, 05}, Vec3{05, 15, 05}, Vec3{-05, 15, 05}},
                // Back face
                {Vec3{05, 0.5, -05}, Vec3{-05, 0.5, -05}, Vec3{-05, 15, -05}},
                {Vec3{05, 0.5, -05}, Vec3{-05, 15, -05}, Vec3{05, 15, -05}},
                // Left face
                {Vec3{-05, 0.5, -05}, Vec3{-05, 0.5, 05}, Vec3{-05, 15, 05}},
                {Vec3{-05, 0.5, -05}, Vec3{-05, 15, 05}, Vec3{-05, 15, -05}},
                // Right face
                {Vec3{05, 0.5, 05}, Vec3{05, 0.5, -05}, Vec3{05, 15, -05}},
                {Vec3{05, 0.5, 05}, Vec3{05, 15, -05}, Vec3{05, 15, 05}},
                // Top face
                {Vec3{-05, 15, 05}, Vec3{05, 15, 05}, Vec3{05, 15, -05}},
                {Vec3{-05, 15, 05}, Vec3{05, 15, -05}, Vec3{-05, 15, -05}},
                // Bottom face
                {Vec3{-05, 0.5, -05}, Vec3{05, 0.5, -05}, Vec3{05, 0.5, 05}},
                {Vec3{-05, 0.5, -05}, Vec3{05, 0.5, 05}, Vec3{-05, 0.5, 05}},
            },color=raylib.GRAY
        },
        Body{
            mesh=[]Triangle{

                // CUBE 5 - Back cube (offset by -1.0 in Z)
                // Front face
                {Vec3{-0.5, -0.5, -0.75}, Vec3{0.5, -0.5, -0.75}, Vec3{0.5, 0.5, -0.75}},
                {Vec3{-0.5, -0.5, -0.75}, Vec3{0.5, 0.5, -0.75}, Vec3{-0.5, 0.5, -0.75}},
                // Back face
                {Vec3{0.5, -0.5, -1.5}, Vec3{-0.5, -0.5, -1.5}, Vec3{-0.5, 0.5, -1.5}},
                {Vec3{0.5, -0.5, -1.5}, Vec3{-0.5, 0.5, -1.5}, Vec3{0.5, 0.5, -1.5}},
                // Left face
                {Vec3{-0.5, -0.5, -1.5}, Vec3{-0.5, -0.5, -0.75}, Vec3{-0.5, 0.5, -0.75}},
                {Vec3{-0.5, -0.5, -1.5}, Vec3{-0.5, 0.5, -0.75}, Vec3{-0.5, 0.5, -1.5}},
                // Right face
                {Vec3{0.5, -0.5, -0.75}, Vec3{0.5, -0.5, -1.5}, Vec3{0.5, 0.5, -1.5}},
                {Vec3{0.5, -0.5, -0.75}, Vec3{0.5, 0.5, -1.5}, Vec3{0.5, 0.5, -0.75}},
                // Top face
                {Vec3{-0.5, 0.5, -0.75}, Vec3{0.5, 0.5, -0.75}, Vec3{0.5, 0.5, -1.5}},
                {Vec3{-0.5, 0.5, -0.75}, Vec3{0.5, 0.5, -1.5}, Vec3{-0.5, 0.5, -1.5}},
                // Bottom face
                {Vec3{-0.5, -0.5, -1.5}, Vec3{0.5, -0.5, -1.5}, Vec3{0.5, -0.5, -0.75}},
                {Vec3{-0.5, -0.5, -1.5}, Vec3{0.5, -0.5, -0.75}, Vec3{-0.5, -0.5, -0.75}},
            },color=raylib.YELLOW
        },*/
    })
    raylib.DisableCursor()

    k :f32 = 0
    for (!raylib.WindowShouldClose()) {
        dt := raylib.GetFrameTime();
        if raylib.IsKeyDown(raylib.KeyboardKey.LEFT ) do c.pos -= Vec3{c.c_right.x,0, c.c_right.z}*dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.RIGHT) do c.pos += Vec3{c.c_right.x,0, c.c_right.z}*dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.A) do c.pos -= Vec3{c.c_right.x,0, c.c_right.z}*dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.D) do c.pos += Vec3{c.c_right.x,0, c.c_right.z}*dt;

        if raylib.IsKeyDown(raylib.KeyboardKey.UP)   do c.pos += Vec3{c.c_forward.x,0, c.c_forward.z}*dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.DOWN) do c.pos -= Vec3{c.c_forward.x,0, c.c_forward.z}*dt;

        if raylib.IsKeyDown(raylib.KeyboardKey.W) do c.pos += Vec3{c.c_forward.x,0, c.c_forward.z}*dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.S) do c.pos -= Vec3{c.c_forward.x,0, c.c_forward.z}*dt;

        if raylib.IsKeyDown(raylib.KeyboardKey.L) do c.c_angle_h += 0.01*dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.R) do c.c_angle_h -= 0.01*dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.H) do c.c_angle_h -= 0.01*dt;

        if raylib.IsKeyDown(raylib.KeyboardKey.SPACE) do c.pos.y -=  dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.LEFT_SHIFT) do c.pos.y +=  dt;

        if raylib.IsKeyDown(raylib.KeyboardKey.U) do c.focal_length += 0.1;
        if raylib.IsKeyDown(raylib.KeyboardKey.P) do c.focal_length -= 0.1;

        mouse_delta := raylib.GetMouseDelta()
        if mouse_delta.x != 0 do c.c_angle_h -= mouse_delta.x/WIDTH * 2*3.1415693;
        if mouse_delta.y != 0 do c.c_angle_v -= mouse_delta.y/HEIGHT * 2*3.1415693;


        update_camera(&c);
        tl := camera_top_left(&c);
        log(c.pos);
        log(c.c_forward);
        log(c.c_angle_h);
        log(tl);
        logln();
        bodies = init_bodies_with_bbox(bodies)
        raylib.BeginDrawing();
        raylib.ClearBackground(raylib.BLACK);

        for body in bodies {
            // fmt.print(body.color);
            texture := render_mesh3(&c, body, body.color);
            // raylib.DrawTexture(texture, 0,0, raylib.WHITE);
            draw_mesh(body.mesh)
        }
        draw_camera(c)
        draw_screen_calcs(c);
        fmt.printf("(%.2f, %.2f, %.2f) fl: %.2f; fps:%.2f\n", c.s_center.x, c.s_center.y, c.s_center.z, c.focal_length, 1/dt)
        fmt.println("end render")
        raylib.EndDrawing();
    }
    raylib.EnableCursor()
}

/*
draw_vec3_text :: proc(pos: raylib.Vector2, v: raylib.Vector3) {
    // Format the vector as a string, e.g. "(x, y, z)"
    text := fmt.tprintf("(%.2f, %.2f, %.2f)", v.x, v.y, v.z)
    defer delete(text)
    c_text := strings.to_cstring(text)
    
    // Draw the text at the given screen position, with font size 20 and white color
    raylib.DrawText(c_text, cast(i32)pos.x, cast(i32)pos.y, 20, raylib.WHITE)
}*/

draw_camera :: proc(c:Camera) {
    scale :f32= SCALE/2
    x := i32(c.pos.x*scale+200)
    y := i32(c.pos.z*scale+200)

    raylib.DrawCircle(x, y, 5, raylib.BLUE);
    dx :i32= x + i32(c.c_forward.x*f32(scale))
    dy :i32= (y + i32(c.c_forward.z*f32(scale)))
    raylib.DrawLine(x,y,dx,dy,raylib.RED);
}

draw_mesh :: proc(mesh: []Triangle) {
    scale :f32= SCALE/2
    for t in mesh {
        for v in t {
            x := i32(v.x*scale+200)
            y := i32(v.z*scale+200)
            raylib.DrawCircle(x, y, 1, raylib.RED);
        }
    }

}
render_mesh :: proc(c:^Camera, mesh: []Triangle) {
    tl := camera_top_left(c)
    k :f32= 1/SCALE;
    fin := false
    for i :i32= 0; i < WIDTH; i+= 1 {
        for j :i32= 0; j < HEIGHT; j+=1 {
            r := Ray{
                c.pos,
                tl + Vec3{f32(i)/SCALE, f32(j)/SCALE,0},
            }

            closest_t := math.inf_f32(1)
            hit_something := false
            for tr in mesh {
                if hit, t := MollerTrumbore(tr[0],tr[1],tr[2], r); hit && t < closest_t {
                    closest_t = t
                    hit_something = true
                    raylib.DrawPixel(i,j,raylib.GRAY);
                }
            }
        }
    }
}
render_mesh2 :: proc(c:^Camera, body: Body, color: Color) {
    tl := camera_top_left(c)

    for i :i32= 0; i < WIDTH; i+= 1 {
        for j :i32= 0; j < HEIGHT; j+=1 {
            // Calculate the target point on the screen plane
            u := f32(i) / f32(WIDTH)   // 0 to 1
            v := f32(j) / f32(HEIGHT)  // 0 to 1
            target := tl + c.s_horizontal * u + c.s_vertical * v;
            direction := normalize3(target - c.pos)
            r := Ray{
                c.pos,
                direction,
            }
            // intersect := ray_box_intersection_simple(r, body.bbox);
            // if !intersect {
            //     // return
            // }

            closest_t := math.inf_f32(1)
            hit_something := false
            for tr in body.mesh {
                if hit, t := MollerTrumbore(tr[0],tr[1],tr[2], r); hit && t < closest_t && t > 0 {
                    closest_t = t
                    hit_something = true
                }
            }

            if hit_something {
                // fmt.print(color);
                // raylib.DrawPixel(i, j, raylib.GRAY)
                raylib.DrawPixel(i, j, color)
            }
        }
    }
}
render_mesh3 :: proc(c:^Camera, body: Body, color: Color) -> raylib.Texture {
    tl := camera_top_left(c)
    slices : [4][]Triangle
    count := len(body.mesh) / 4
    remainder := len(body.mesh) % 4

    // Create image buffer
    image := raylib.GenImageColor(WIDTH, HEIGHT, raylib.BLACK)
    image.format = raylib.PixelFormat.UNCOMPRESSED_R8G8B8
    defer raylib.UnloadImage(image)
    
    ThreadData :: struct {
        mesh: []Triangle, 
        color: Color, 
        i: i32,
        tl: Vec3, 
        c: Camera,
        image_data: [^]Color, // Raw pointer to image data
    }

    threads := [dynamic]^thread.Thread{}
    defer delete(threads)

    for i :i32= 0; i < WIDTH; i+= 1 {
        f:: proc(data: rawptr) {
            tdata := cast(^ThreadData)data
            for j :i32= 0; j < HEIGHT; j+=1 {
                // Calculate the target point on the screen plane
                u := f32(tdata.i) / f32(WIDTH)   // 0 to 1
                v := f32(j) / f32(HEIGHT)  // 0 to 1
                target := tdata.tl + tdata.c.s_horizontal * u + tdata.c.s_vertical * v;
                direction := normalize3(target - tdata.c.pos)
                r := Ray{
                    tdata.c.pos,
                    direction,
                }
                closest_t := math.inf_f32(1)
                hit_something := false
                for tr in tdata.mesh {
                    if hit, t := MollerTrumbore(tr[0],tr[1],tr[2], r); hit && t < closest_t && t > 0 {
                        closest_t = t
                        hit_something = true
                    }
                }
                if hit_something {
                    // fmt.println("Hit");
                    index := j * WIDTH + tdata.i
                    tdata.image_data[index] = tdata.color
                    // fmt.println(index, tdata.color);
                    // raylib.DrawPixel(tdata.i, j, tdata.color)
                }
            }
            free(data)
        }
        data := new(ThreadData)
        data^ = ThreadData{
            mesh=body.mesh,
            color=body.color,
            i=i,
            tl=tl,
            c=c^,
            image_data = cast([^]Color)image.data
        }
        t := thread.create_and_start_with_data(data, f)
        append(&threads, t)
    }
    // fmt.println(len(threads));
    for t in threads {
        thread.join(t)
    }
    // Destroy threads
    for t in threads {
        thread.destroy(t)
    }
        // Get raw access to image data
    image_data := cast([^]Color)image.data
    for i in 0..<WIDTH {
        for j in 0..<HEIGHT {
            raylib.DrawPixel(i32(i), i32(j), image_data[j*WIDTH + i]);
        }

    }
    return raylib.Texture2D{}
    /*
    // Check if any pixels were actually set
    pixels_set := 0
    background := raylib.BLACK
    for i in 0..<(WIDTH * HEIGHT) {
        pixel := image_data[i]
        if pixel.r != background.r || 
           pixel.g != background.g || 
           pixel.b != background.b || 
           pixel.a != background.a {
            pixels_set += 1
        }
    }
    
    fmt.printf("Total pixels set: %d out of %d\n", pixels_set, WIDTH * HEIGHT)
    
    if pixels_set == 0 {
        fmt.println("WARNING: No pixels were set by raytracing!")
    }
    return raylib.Texture2D{}
    */
}

check_uniform_texture :: proc(texture: raylib.Texture2D) -> bool {
    img := raylib.LoadImageFromTexture(texture)
    defer raylib.UnloadImage(img)

    pixels := raylib.LoadImageColors(img)
    defer raylib.UnloadImageColors(pixels)
    first := pixels[0]
    total := img.width * img.height
    width  := img.width
    height := img.height


    totalc := raylib.Color{0, 0, 0,0}
    pixel_count := width * height

    for i in 0..<pixel_count {
        totalc.r += pixels[i].r
        totalc.g += pixels[i].g
        totalc.b += pixels[i].b
        totalc.a += pixels[i].a
    }

    avg := raylib.Color{
        u8(totalc.r / u8(pixel_count)),
        u8(totalc.g / u8(pixel_count)),
        u8(totalc.b / u8(pixel_count)),
        u8(totalc.a / u8(pixel_count)),
    }

    fmt.println("IN CHECK UNIFORM:", avg);


    if first.r == 0 && first.b == 0 && first.g == 0 {
        fmt.println("BLANK!!!!! ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK ITS BLANK")
    }
    for i in 1..<total {
        if pixels[i] != first {
            return false // Not uniform
        }
    }


    return true // All pixels are the same
}


draw_screen_calcs  :: proc(c:Camera) {
    cam := c
    scale := SCALE/3
    x0 :i32= WIDTH / 2
    y0 :i32= HEIGHT / 2
    m := []Vec3{
        camera_top_left(&cam),
        camera_top_right(&cam),
        camera_bottom_left(&cam),
        camera_bottom_right(&cam),
    }
    for i in m {
        log(i)
        raylib.DrawCircle(x0+i32(i.x*SCALE), y0+i32(i.y*SCALE), 2,raylib.WHITE);
    }
    logln()
    raylib.DrawCircle(x0+i32(c.s_center.x*SCALE), y0+i32(c.s_center.y*SCALE), 2,raylib.WHITE);

}
average_texture_color :: proc(texture: raylib.Texture2D) -> Color {
    image := raylib.LoadImageFromTexture(texture)
    defer raylib.UnloadImage(image)

    width  := image.width
    height := image.height

    pixels := raylib.LoadImageColors(image)
    defer raylib.UnloadImageColors(pixels)

    total := raylib.Color{0, 0, 0,0}
    pixel_count := width * height

    for i in 0..<pixel_count {
        total.r += pixels[i].r
        total.g += pixels[i].g
        total.b += pixels[i].b
        total.a += pixels[i].a
    }

    avg := raylib.Color{
        u8(total.r / u8(pixel_count)),
        u8(total.g / u8(pixel_count)),
        u8(total.b / u8(pixel_count)),
        u8(total.a / u8(pixel_count)),
    }

    return avg
}

