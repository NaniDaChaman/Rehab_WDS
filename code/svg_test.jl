using NativeSVG
Drawing() do
    g(stroke_linecap="butt", stroke_miterlimit="4", stroke_width="3.0703125") do
        circle(cx="20", cy="20", r="16", stroke="#cb3c33", fill="#d5635c")
        circle(cx="40", cy="56", r="16", stroke="#389826", fill="#60ad51")
        circle(cx="60", cy="20", r="16", stroke="#9558b2", fill="#aa79c1")
    end
end