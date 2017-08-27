function [a b c] = switch_master(T_sl, T_sg, T_boil, m_sl, m_s)
	global state;
	m_buffer = 0.05*m_s;
	if (state(1) == 0 && state(2) == 1 && state(3) == 0)
		if (T_sl > T_boil)
			state = [1, 0, 1];
		end
	elseif (state(1) == 1 && state(2) == 0 && state(3) == 1)
		if (m_sl < m_buffer)
			state = [0, 0, 1];
		elseif (m_sl > (m_s - m_buffer))
			state = [0, 1, 0];
		end
	elseif (state(1) == 0 && state(2) == 0 && state(3) == 1)
		if (T_sg < T_boil)
			state = [1, 0, 1];
		end
	end
	a = state(1);
	b = state(2);
	c = state(3);
end